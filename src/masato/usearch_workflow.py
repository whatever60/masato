#!/usr/bin/env python
import argparse
import glob
import io
import os
import re
import subprocess
import shutil
import gzip
from typing import IO
from io import StringIO
import tempfile
from concurrent.futures import as_completed
import json
from collections import defaultdict

import numpy as np
from scipy import sparse as ss
import pandas as pd
from loky import get_reusable_executor
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import biom
import anndata as ad
from tqdm.auto import tqdm

from masato.utils import read_table, write_table


def merge_pairs(
    fastq_dir: str, fastq_out: str, num_threads: int, backend: str = "vsearch"
) -> None:
    # merge paired end reads
    output_dir = os.path.dirname(fastq_out)
    os.makedirs(output_dir, exist_ok=True)
    if backend == "usearch":
        subprocess.run(
            "usearch11 -fastq_mergepairs "
            f"{fastq_dir}/*_R1.fastq "
            "-relabel @ "
            f"-fastqout {output_dir}/merged.fastq "
            "-fastq_maxdiffs 25 "
            "-fastq_minovlen 50 "
            f"-report {output_dir}/merge_report.txt "
            f"-tabbedout {output_dir}/merge.tsv",
            shell=True,
        )
    elif backend == "vsearch":
        # vsearch cannot merge multiple files or perform relabeling at once, so we first
        # relabel with seqtk and then merge samples, and finally merge pairs with
        # vsearch.
        # Besides, vsearch does not have tsv or report output
        os.makedirs(f"{output_dir}/temp", exist_ok=True)
        fs = []
        ps = []

        # `fastq_dir` could be a path regex so that we can easily do merging for a
        # subset of samples in a directory.
        if os.path.isdir(fastq_dir):
            files = sorted(glob.glob(f"{fastq_dir}/*.fastq"))
            if not files:
                raise ValueError(f"No fastq files found in {fastq_dir}")
        else:
            files = [i for i in sorted(glob.glob(fastq_dir)) if i.endswith(".fastq")]
            if not files:
                raise ValueError(f"No fastq files found matching {fastq_dir}")
        if not all(f.endswith("_R1.fastq") or f.endswith("_R2.fastq") for f in files):
            raise ValueError(
                "Make sure to rename all fastq files of interest to end with _R1.fastq or _R2.fastq"
            )
        non_pair_files = []
        for f in files:
            if f.endswith("_R1.fastq"):
                if f"{f[:-9]}_R2.fastq" not in files:
                    non_pair_files.append(f)
            elif f.endswith("_R2.fastq"):
                if f"{f[:-9]}_R1.fastq" not in files:
                    non_pair_files.append(f)
        if non_pair_files:
            raise ValueError(
                f"Make sure that each sample has both _R1.fastq and _R2.fastq files. "
                f"The following files are not paired: {non_pair_files}"
            )

        batch = 1000
        for i in range(0, len(files), batch):
            file_list = files[i : i + batch]
            for f in file_list:
                f_base = os.path.basename(f)
                sample_name = f_base.split("_")[0]
                # use Popen for non-blocking execution
                outfile = open(f"{output_dir}/temp/{f_base}", "w")
                process = subprocess.Popen(
                    ["seqtk", "rename", f, f"sample={sample_name} "],
                    stdout=outfile,
                )
                fs.append(outfile)
                ps.append(process)
            [p.wait() for p in tqdm(ps)]
            [f.close() for f in fs]

        subprocess.run(
            f"cat {output_dir}/temp/*_R1.fastq > {output_dir}/temp_R1.fastq", shell=True
        )
        subprocess.run(
            f"cat {output_dir}/temp/*_R2.fastq > {output_dir}/temp_R2.fastq", shell=True
        )
        subprocess.run(
            [
                "vsearch",
                "--fastq_mergepairs",
                f"{output_dir}/temp_R1.fastq",
                "--reverse",
                f"{output_dir}/temp_R2.fastq",
                "--fastqout",
                fastq_out,
                # "--fastaout",
                # f"{output_dir}/merged.fa",
                "--fastq_maxdiffs",
                "25",
                "--fastq_minovlen",
                "50",
                "--threads",
                str(num_threads),
            ]
        )
        os.remove(f"{output_dir}/temp_R1.fastq")
        os.remove(f"{output_dir}/temp_R2.fastq")
        shutil.rmtree(f"{output_dir}/temp")


def qc(
    input_fastq: str,
    output_fasta: str | None = None,
    num_threads: int = 16,
    backend: str = "vsearch",
) -> str:
    """Output fasta is default to be modified from input fastq, but can be specified.
    For example, if input fastq file is /data/merged.fq.gz, then output fasta will be
    /data/merged.filtered.fq.gz, unless otherwise specified.
    """
    if output_fasta is None:
        output_dir = os.path.dirname(input_fastq)
        filename = os.path.basename(input_fastq).split(".")[0]
        output_fastx = os.path.join(output_dir, f"{filename}.filtered.fq")
    else:
        output_fastx = output_fasta.rstrip(".gz")

    if backend == "usearch":
        num_splits = None
        # Split the merged FASTQ file into multiple pieces
        total_lines = int(subprocess.check_output(["wc", "-l", input_fastq]).split()[0])
        lines_per_file = ((total_lines + num_splits - 1) // num_splits // 4) * 4
        base_dir = os.path.dirname(input_fastq)
        if num_splits > 1:
            subprocess.run(
                [
                    "split",
                    "-l",
                    str(lines_per_file),
                    "-d",
                    "--additional-suffix=.fastq",
                    input_fastq,
                    os.path.join(base_dir, "_split_"),
                ]
            )

            # Perform QC on each split FASTQ file
            with open(output_fastx, "w") as outfile:
                for split_file in glob.glob(os.path.join(base_dir, "_split_*.fastq")):
                    subprocess.run(
                        [
                            "usearch11",
                            "-fastq_filter",
                            split_file,
                            "-fastaout",
                            "/dev/stdout",
                            "-fastq_maxee",
                            "0.5",
                            "-fastq_minlen",
                            "200",
                            "-fastq_maxns",
                            "0",
                            "-relabel",
                            "filtered",
                            "-threads",
                            str(num_threads),
                        ],
                        stdout=outfile,
                    )
                    os.remove(split_file)
        else:
            subprocess.run(
                [
                    "usearch11",
                    "-fastq_filter",
                    input_fastq,
                    "-fastaout",
                    output_fastx,
                    "-fastq_maxee",
                    "0.5",
                    "-fastq_minlen",
                    "200",
                    "-fastq_maxns",
                    "0",
                    "-relabel",
                    "filtered",
                    "-threads",
                    str(num_threads),
                ]
            )
    elif backend == "vsearch":
        # we don't need to split anymore since there is no memory limit.
        subprocess.run(
            [
                "vsearch",
                "--fastq_filter",
                input_fastq,
                "--fastqout",
                output_fastx,
                "--fastq_maxee",
                "0.5",
                "--fastq_minlen",
                "100",
                "--fastq_maxns",
                "0",
                "--relabel",
                "filtered",
                # "--threads",
                # str(num_threads),
            ]
        )
        subprocess.run(["pigz", "-f", output_fastx])
        output_fastx += ".gz"
    return output_fastx


def split_fastq(args):
    # Calculate the number of lines in the FASTQ file (each read consists of 4 lines)
    total_lines = int(subprocess.check_output(["wc", "-l", args.input_file]).split()[0])
    # Calculate the number of lines for each split file, ensuring that the total number of lines is divisible by 4
    lines_per_file = ((total_lines + args.num_splits - 1) // args.num_splits // 4) * 4

    # Split the file using the `split` command
    subprocess.run(
        [
            "split",
            "-l",
            str(lines_per_file),
            "-d",
            "--additional-suffix=.fastq",
            args.input_file,
            os.path.join(args.output_dir, "split_"),
        ]
    )

    print(f"FASTQ file has been split into {args.num_splits} pieces.")


def add_depth_to_metadata(fastq_dir: str, metadata_path: str) -> None:
    # Generate metadata for the input FASTQ files

    # data = []
    # for f in os.listdir(fastq_dir):
    #     if f.endswith("_R1.fastq"):
    #         sample_name = f.replace("_R1.fastq", "")
    #         count = sum(1 for line in open(os.path.join(fastq_dir, f))) / 4
    #         data.append([sample_name, count])
    command = f"seqkit stats {fastq_dir}/*_R1.fastq -T"
    result = subprocess.run(command, capture_output=True, text=True, shell=True)

    # Check if the command was successful
    if result.returncode != 0:
        print("Error executing command:", result.stderr)
        exit()

    # Split the captured output into lines
    lines = result.stdout.splitlines()

    # Use pandas to read the tab-separated output into a DataFrame
    df_rc = pd.read_csv(io.StringIO("\n".join(lines)), sep="\t")
    df_rc["sample"] = df_rc.file.apply(
        lambda x: os.path.basename(x).replace("_R1.fastq", "")
    )
    df_rc = df_rc[["sample", "num_seqs"]].set_index("sample")
    df_rc.columns = ["read_count"]

    df = pd.read_table(metadata_path, index_col="sample")
    if "read_count" in df.columns:
        df.drop("read_count", axis=1, inplace=True)
    # raise error if there are samples in the metadata that are not in the fastq
    # directory, and raise warning if there are samples in the fastq directory that
    # are not in the metadata (these samples will be ignored)
    extra_samples = df_rc.index.difference(df.index)
    if len(extra_samples) > 0:
        print(
            f"WARNING: These samples in the fastq directory are not present in the metadata: {extra_samples}"
        )
    extra_samples = df.index.difference(df_rc.index)
    if len(extra_samples) > 0:
        raise ValueError(
            f"These samples in the metadata are not present in the fastq directory: {extra_samples}"
        )
    else:
        print("Great! All samples in the metadata are present in the fastq directory.")
    df_res = pd.merge(df, df_rc, left_index=True, right_index=True, how="inner")
    df_res.to_csv(metadata_path, sep="\t")


def subsample(
    fastq_dir: str,
    metadata_path: str,
    output_dir: str,
    output_metadata_path: str,
    num_subsamples: int,
    mode: str,
    seed: int = 100,
    metadata_only: bool = False,
) -> None:
    # Load the metadata CSV file into a pandas DataFrame
    df_meta = pd.read_table(metadata_path, index_col="sample")

    if mode == "reference":
        df_meta["read_count_subsample"] = df_meta["reference"].map(
            df_meta["read_count"]
        )
    elif mode == "min":
        df_meta["read_count_subsample"] = df_meta.read_count.min() - 1
    else:
        raise ValueError(f"Invalid mode: {mode}")

    os.makedirs(output_dir, exist_ok=True)

    # Begin subsampling
    subsample_metadata = []
    for sample_name, row in tqdm(df_meta.iterrows(), total=len(df_meta)):
        read_count = min(int(row["read_count_subsample"]), int(row["read_count"]))
        num_subs = (
            num_subsamples
            if int(row["read_count_subsample"]) < int(row["read_count"])
            else 1
        )
        del row["read_count_subsample"], row["read_count"]
        ps = []
        for iteration in range(1, num_subs + 1):
            # Subsample R1 and R1 and output to the subsampled directory
            subsample_metadata.append(
                {
                    "sample": f"{sample_name}-subsampled-{iteration}",
                    "original_sample": sample_name,
                    "read_count": read_count,
                }
                | row.to_dict()
            )
            if metadata_only:
                continue
            # Increment the seed value
            seed += 1
            for suffix in ["_R1.fastq", "_R2.fastq"]:
                fastq_file = os.path.join(fastq_dir, f"{sample_name}{suffix}")
                if not os.path.isfile(fastq_file):
                    print(f"WARNING: {fastq_file} does not exist. Skipping.")
                    continue
                # can also use seqtk, but will need to redirect stdout to file
                p = subprocess.Popen(  # redirect stdout and stderr to /dev/null
                    [
                        "vsearch",
                        "--fastx_subsample",
                        fastq_file,
                        "--sample_size",
                        str(read_count),
                        "-randseed",
                        str(seed),
                        "--fastqout",
                        os.path.join(
                            output_dir, f"{sample_name}-subsampled-{iteration}{suffix}"
                        ),
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
                ps.append(p)
        [p.wait() for p in ps]
    pd.DataFrame(subsample_metadata).to_csv(output_metadata_path, sep="\t", index=False)


def db_construct(
    input_fastx: str,
    output_fasta: str | None,
    num_threads: int,
    backend: str = "vsearch",
) -> str | None:
    """Output fasta is default to be modified from input fastq, but can be specified.
    For example, if input fastq file is /data/merged.filtered.fq.gz, then output fasta
    will be /data/merged.uniq.fa.gz, unless otherwise specified.
    """
    if output_fasta is None:
        output_dir, base_name = os.path.split(input_fastx)
        base_name_noext = base_name.split(".")[0]
        output_path = f"{output_dir}/{base_name_noext}.uniq.fa"
    else:
        output_path = output_fasta.rstrip(".gz")
        output_dir = os.path.dirname(output_path)
        base_name_noext = os.path.basename(output_path).split(".")[0]

    if backend == "usearch":
        subprocess.run(
            [
                "usearch11",
                "-fastx_uniques",
                input_fastx,
                "-fastaout",
                f"{output_dir}/{base_name_noext}.uniq.fa",
                "-sizeout",
            ],
            check=True,
        )
        subprocess.run(
            [
                "usearch11",
                "-sortbysize",
                f"{output_dir}/{base_name_noext}.uniq.fa",
                "-fastaout",
                f"{output_dir}/{base_name_noext}.sorted.fa",
                "-minsize",
                "2",
            ],
            check=True,
        )
        subprocess.run(
            [
                "usearch11",
                "-search_exact",
                input_fastx,
                "-db",
                f"{output_dir}/{base_name_noext}.sorted.fa",
                "-strand",
                "both",
                "-notmatched",
                f"{output_dir}/{base_name_noext}.unmatched.fa",
                "-t",
                str(num_threads),
            ],
            check=True,
        )
    elif backend == "vsearch":
        subprocess.run(
            [
                "vsearch",
                "--fastx_uniques",
                input_fastx,
                "--fastaout",
                output_path,
                "--sizeout",
            ]
        )
        subprocess.run(["pigz", "-f", output_path])
        output_path += ".gz"
        return output_path


def cluster_uparse(
    input_fastq: str, db_fasta: str, output_dir: str, num_threads: int
) -> None:
    if output_dir is None:
        output_dir = os.path.dirname(input_fastq)
    os.makedirs(output_dir, exist_ok=True)
    subprocess.run(
        [
            "usearch11",
            "-cluster_otus",
            db_fasta,
            "-otus",
            f"{output_dir}/uparse_otu.fa",
            "-relabel",
            "OTU",
        ],
        check=True,
    )
    subprocess.run(
        [
            "usearch10",
            "-usearch_global",
            input_fastq,
            "-db",
            f"{output_dir}/uparse_otu.fa",
            "-strand",
            "both",
            "-id",
            "0.97",
            "-otutabout",
            f"{output_dir}/uparse_otu.tsv",
            "-biomout",
            f"{output_dir}/uparse_otu.biom",
            "-threads",
            str(num_threads),
            "-sample_delim",
            ".",
        ],
        check=True,
    )
    subprocess.run(
        [
            "usearch11",
            "-calc_distmx",
            f"{output_dir}/uparse_otu.fa",
            "-tabbedout",
            f"{output_dir}/uparse_otu.dist",
            "-maxdist",
            "0.2",
            "-termdist",
            "0.3",
            "-threads",
            str(num_threads),
        ],
        check=True,
    )
    subprocess.run(
        [
            "usearch11",
            "-cluster_aggd",
            f"{output_dir}/uparse_otu.dist",
            "-treeout",
            f"{output_dir}/uparse_otu.tree",
            "-clusterout",
            f"{output_dir}/uparse_otu.clust",
            "-id",
            "0.8",
            "-linkage",
            "min",
        ],
        check=True,
    )


def unoise3(
    input_fastq: str | IO[str] | gzip.GzipFile,
    output_fasta: str | IO[str] | gzip.GzipFile | None,
    min_size: int,
    alpha: float = 2.0,
    min_length: int = 0,
    relabel_prefix: str | None = None,
    size_out: bool = False,
    num_threads: int = 16,
    stderr: bool = True,
) -> None | str:
    stderr_ = subprocess.DEVNULL if not stderr else None
    input_arg, stdin = _decide_io_arg(input_fastq)
    output_arg, stdout = _decide_io_arg(output_fasta)

    args_qc = [
        "vsearch",
        "--fastq_filter",
        input_arg,
        "--fastqout",
        "-",
        "--fastq_maxee",
        "0.5",
        "--fastq_minlen",
        str(min_length),
        "--fastq_maxns",
        "0",
        "--relabel",
        "filtered",
        # "--threads",
        # str(num_threads),
    ]
    # print_command(args_qc)
    qc_proc = subprocess.Popen(
        args_qc,
        stdin=stdin,
        stdout=subprocess.PIPE,
        stderr=stderr_,
    )

    args_uniq = [
        "vsearch",
        "--fastx_uniques",
        "-",
        "--fastaout",
        "-",
        "--sizeout",
    ]
    # print_command(args_uniq)
    uniq_proc = subprocess.Popen(
        args_uniq, stdin=qc_proc.stdout, stdout=subprocess.PIPE, stderr=stderr_
    )

    args_unoise3 = [
        "vsearch",
        "--cluster_unoise",
        "-",
        "--minsize",
        str(min_size),
        "--unoise_alpha",
        str(alpha),
        "--sizeout",
        "--centroids",
        "-",
        "--threads",
        str(num_threads),
    ]
    # print_command(args_unoise3)
    unoise3_proc = subprocess.Popen(
        args_unoise3, stdin=uniq_proc.stdout, stdout=subprocess.PIPE, stderr=stderr_
    )

    args_uchime3 = [
        "vsearch",
        "--uchime3_denovo",
        "-",
        "--nonchimeras",
        output_arg,
        "--relabel",
        relabel_prefix or "ZOTU",
    ]
    if size_out:
        args_uchime3.append("--sizeout")
    uchime3_proc = subprocess.Popen(
        args_uchime3, stdin=unoise3_proc.stdout, stdout=stdout, stderr=stderr_
    )

    if qc_proc.stdout is not None:
        qc_proc.stdout.close()
    if uniq_proc.stdout is not None:
        uniq_proc.stdout.close()
    if unoise3_proc.stdout is not None:
        unoise3_proc.stdout.close()
    if stdin is subprocess.PIPE:
        if isinstance(input_fastq, str):
            qc_proc.communicate(input_fastq.encode())
        else:
            # If input_fastq is a file-like object, we shouldn't reach here
            # because _decide_io_arg would have set stdin to input_fastq, not PIPE
            raise TypeError(
                f"Expected input_fastq to be str when stdin is PIPE, got {type(input_fastq)}"
            )
    uchime3_out, _ = uchime3_proc.communicate()
    if stdout is subprocess.PIPE:
        return uchime3_out.decode()


def search_global(
    input_fastq: str | gzip.GzipFile | IO[str] | None,
    zotu_fasta: str | gzip.GzipFile | IO[str] | None,
    output_tsv: str | gzip.GzipFile | IO[str] | None,
    not_matched_fasta: str | gzip.GzipFile | IO[str] | None = None,
    id_: float = 0.97,
    num_threads: int = 8,
    stderr: bool = True,
) -> None | str:
    """
    When `output_csv` is None, output is in captured in the stdout of returned process.
    When `not_matched_fasta` is None, unmatched query sequences are not saved as fasta.
    """
    fastq_arg, stdin = _decide_io_arg(input_fastq)
    db_arg, stdin_db = _decide_io_arg(zotu_fasta)
    tsv_arg, stdout = _decide_io_arg(output_tsv)
    if not_matched_fasta is None:
        not_matched_arg, stdout_not_matched = None, None
    else:
        not_matched_arg, stdout_not_matched = _decide_io_arg(not_matched_fasta)
    if stdin is not None and stdin_db is not None:
        raise ValueError(
            "`input_fastq` and `zotu_fasta` cannot both be file-like objects."
        )
    if stdout is not None and stdout_not_matched is not None:
        raise ValueError(
            "`output_tsv` and `not_matched_fasta` cannot both be file-like objects."
        )
    stderr = subprocess.DEVNULL if not stderr else None
    # https://github.com/torognes/vsearch/issues/552
    # https://github.com/torognes/vsearch/issues/467
    # https://github.com/torognes/vsearch/issues/392
    args_search = [
        "vsearch",
        "--usearch_global",
        fastq_arg,
        "--db",
        db_arg,
        "--strand",
        "plus",
        "--id",
        str(id_),
        "--maxaccepts",
        "10",
        "--maxhits",
        "1",
        "--notrunclabels",
        "--no_progress",
        "--threads",
        str(num_threads),
    ]
    if tsv_arg == "-":
        args_search.extend(["--otutabout", "-"])
    elif tsv_arg.endswith(".tsv") or tsv_arg.endswith(".tsv.gz"):
        args_search.extend(["--otutabout", tsv_arg])
    elif output_tsv.endswith(".biom"):
        args_search.extend(["--biomout", tsv_arg])
    else:
        raise ValueError(
            f"Invalid output file extension {output_tsv}. Must be one of .tsv, .biom"
        )
    if not_matched_arg is not None:
        args_search.extend(["--notmatched", not_matched_arg])
    sg_proc = subprocess.Popen(
        args_search,
        stdin=stdin,
        stdout=stdout,
        stderr=stderr,
    )
    sg_proc_out, sg_proc_err = sg_proc.communicate(
        input_fastq.encode() if stdin is subprocess.PIPE else None
    )
    if stdout is subprocess.PIPE:
        return sg_proc_out.decode()


def csr_vappend(a: ss.csr_matrix, b: ss.csr_matrix) -> ss.csr_matrix:
    """Takes in 2 csr_matrices and appends the second one to the bottom of the first one.
    Much faster than scipy.sparse.vstack but assumes the type to be csr and overwrites
    the first matrix instead of copying it. The data, indices, and indptr still get copied."""

    a.data = np.hstack((a.data, b.data))
    a.indices = np.hstack((a.indices, b.indices))
    a.indptr = np.hstack((a.indptr, (b.indptr + a.nnz)[1:]))
    a._shape = (a.shape[0] + b.shape[0], b.shape[1])
    return a


def search_global_add_unknown(
    input_fastq: str,
    zotu_fasta: str,
    output_path: str,
    id_: float,
    unknown_name: str,
    num_threads: int,
) -> None:
    if not os.path.isfile(zotu_fasta):
        raise ValueError(f"ZOTU fasta file (database) {zotu_fasta} does not exist")
    if not os.path.isfile(input_fastq):
        raise ValueError(f"Query fastq file {input_fastq} does not exist")
    if output_path is None:
        output_dir = os.path.dirname(input_fastq)
        output_tsv = os.path.join(output_dir, "unoise3_zotu.biom")
    else:
        output_dir = os.path.dirname(output_path)
        os.makedirs(output_dir, exist_ok=True)
    output_not_matched = os.path.join(output_dir, "notmatched.fa")

    search_global(
        input_fastq, zotu_fasta, output_path, output_not_matched, id_, num_threads
    )
    _add_unknown(output_path, output_not_matched, unknown_name, output_path)


def get_sample_name_from_header(header: str) -> str:
    """Extract "sample=" from one of the header fields"""
    for field in re.split(r"[;\s]+", header.strip()):
        if field.startswith("sample="):
            return field.split("=")[1]
    else:
        raise ValueError(f"Header {header} does not contain 'sample=' field.")

def _add_unknown(
    input_path: str, output_not_matched: str, unknown_name: str, output_path: str
) -> None:
    # add a ZOTU_UNKNOWN to the zotu table by counting the number of reads that do not
    # matched to anything in database
    zotu_table: ad.AnnData = read_table(input_path, return_type="adata")
    matrix_data: ss.csr_matrix = zotu_table.X
    samples = zotu_table.obs_names
    counter = {i: 0 for i in samples}
    for record in SeqIO.parse(output_not_matched, "fasta"):
        counter[get_sample_name_from_header(record.description)] += 1

    matrix_data = ss.hstack(
        [matrix_data, ss.csr_matrix([counter[i] for i in samples]).T], format="csr"
    )
    obs = zotu_table.obs
    var_ = pd.concat(
        [zotu_table.var, pd.DataFrame(index=[unknown_name], dtype=matrix_data.dtype)]
    )
    zotu_table = ad.AnnData(matrix_data, obs=obs, var=var_)
    write_table(zotu_table, output_path)


def _decide_io_arg(arg) -> tuple:
    if arg is None:  # PIPE as IO
        return "-", subprocess.PIPE
    elif isinstance(arg, str):
        if "\n" in arg:  # string literal corresponding to file content
            return "-", subprocess.PIPE
        else:  # file path
            arg = os.path.abspath(arg)
            if (
                not os.path.isfile(arg)
                and not os.path.isdir(arg)
                and not os.path.isdir(os.path.dirname(arg))
            ):
                raise ValueError(
                    f"{arg} does not appear to be stirng liternal and is also not a file path."
                )
            return arg, None
    elif isinstance(arg, StringIO):
        return "-", subprocess.PIPE
    else:
        if (
            not isinstance(arg, gzip.GzipFile)
            and not isinstance(arg, IO[str])
            and not isinstance(arg, subprocess.Popen)
        ):
            raise ValueError(f"Invalid argument {arg} of type {type(arg)}")
        return "-", arg


def _workflow_one_sample(
    seqs_sample: list[str],
    min_size: int,
    alpha: float,
    prefix: str | None = None,
    search: bool = True,
) -> tuple[list[str], list[str], list[int]]:
    num_qs = len(seqs_sample)
    input_fastq = "".join(seqs_sample)
    db_fasta = tempfile.NamedTemporaryFile()
    unoise3(
        input_fastq,
        db_fasta.name,
        min_size=min_size,
        alpha=alpha,
        relabel_prefix=prefix,
        size_out=False,
        num_threads=1,
        stderr=False,
    )
    try:
        zotus_name, zotus_seq = zip(
            *(
                (record.id, str(record.seq))
                for record in SeqIO.parse(db_fasta.name, "fasta")
            )
        )
    except ValueError:
        zotus_name, zotus_seq = [], []
    else:
        zotus_name, zotus_seq = list(zotus_name), list(zotus_seq)
    if zotus_name and search:
        count_tsv = search_global(
            input_fastq, db_fasta.name, None, num_threads=1, stderr=False
        )
        counts = pd.read_table(StringIO(count_tsv), index_col=0).iloc[:, 0].tolist()
        num_unknown = num_qs - sum(counts)
        if num_unknown > 0:
            zotus_name.append("#UNKNOWN")
            zotus_seq.append("#UNKNOWN")
            counts.append(num_unknown)
    else:
        counts = []
        if num_qs > 0:
            zotus_name = ["#UNKNOWN"]
            zotus_seq = ["#UNKNOWN"]
            counts = [num_qs]
    db_fasta.close()
    return zotus_name, zotus_seq, counts


def workflow_per_sample(
    input_fastq: str,
    output_path: str,
    min_size: int,
    alpha: float,
    prefix: str | None = None,
    num_threads: int = 8,
    search: bool = True,
) -> None:
    """Split input fastq file by sample and run a separate unoise3 workflow for each
        sample, so that each sample has one set of ZOTUs and counts.

    We assume sequences in the input fastq file are named like `@sample=<sample1>` and
        sequences from the same sample form consecutive blocks.

    If output path ends with json, save results in a json file like:
        `{<sample>: {"zotus": ["ATCG", "GCTA", ...], "counts": [10, 2, ...]}, ...}`.
    If output path ends with fasta, save results in a fasta file file like:
        ```
        >ATCG
        ATCG
        ```
        where sequence names are just sequences themselves. Identical sequences appear
        only once.
    """
    from masato.utils import smart_open

    current_sample = None
    fastq_for_current_sample = []
    future2sample = {}
    samples = set()
    # sample2future = {}
    results = {}

    executor = (
        get_reusable_executor(max_workers=num_threads, kill_workers=True, reuse=False)
        if num_threads > 1
        else None
    )
    f: IO[str] = smart_open(input_fastq, "r")
    for entry in zip(f, f, f, f):
        # Parse the header line to extract the sample name
        header: str = entry[0]
        sample_name = get_sample_name_from_header(header)
        if sample_name != current_sample and fastq_for_current_sample:
            if sample_name in samples:
                raise ValueError(
                    "Sequences in fastq file are not grouped by sample, found "
                    f"{sample_name} both before {current_sample} and after"
                )
            if executor:
                future = executor.submit(
                    _workflow_one_sample,
                    fastq_for_current_sample,
                    min_size,
                    alpha,
                    prefix,
                    search,
                )
                future2sample[future] = current_sample
                # sample2future[current_sample] = future
            else:
                results[current_sample] = _workflow_one_sample(
                    fastq_for_current_sample, min_size, alpha, prefix, search
                )
            fastq_for_current_sample = []  # Reset for the next sample
            samples.add(current_sample)
        current_sample = sample_name
        fastq_for_current_sample.append("".join(entry))
    f.close()
    # Don't forget to process the last sample
    if fastq_for_current_sample:
        if executor:
            future = executor.submit(
                _workflow_one_sample,
                fastq_for_current_sample,
                min_size,
                alpha,
                prefix,
                search,
            )
            future2sample[future] = current_sample
            # sample2future[current_sample] = future
            for future in tqdm(as_completed(future2sample), total=len(future2sample)):
                # start the task and save the future object
                sample = future2sample[future]
                results[sample] = future.result()
            results = {s: results[s] for s in samples}
        else:
            results[current_sample] = _workflow_one_sample(
                fastq_for_current_sample, min_size, alpha, prefix, search
            )
    # Assuming each task's result could be aggregated into a final result dictionary
    # This step depends on how you implement the process_queue and task results handling
    # You would collect results from each completed task here
    # This could involve collecting returned values or reading from a shared resource
    ext = os.path.splitext(output_path)[1]
    if ext == ".json":
        results = {s: {"zotus": z, "counts": c} for s, (_, z, c) in results.items()}
        with open(output_path, "w") as outfile:
            json.dump(results, outfile, indent=4)
    elif ext in set([".fa", ".fasta", ".fna"]):
        prefix = prefix or "ZOTU"
        zotus = {}
        for zns, zss, counts in results.values():
            for zn, zs, count in zip(zns, zss, counts):
                if zn == "#UNKNOWN":
                    continue
                zotus[zs] = zotus.get(zs, 0) + count
        # sort by frequency
        records = [
            SeqRecord(Seq(z), id=f"{prefix}{i}", description=str(zotus[z]))
            for i, z in enumerate(sorted(zotus, key=zotus.get, reverse=True), 1)
        ]
        with open(output_path, "w") as outfile:
            SeqIO.write(records, outfile, "fasta")
    else:
        raise ValueError(
            f"Invalid output file extension {ext}. Must be one of .json, .fa, .fasta, .fna"
        )


def aggregate_samples(
    input_json: str, output_fasta: str, output_count: str, prefix: str | None = None
) -> None:
    """Take the json file generated by workflow_per_sample and aggregate the results into
    fasta file of ZOTUs and count file.

    ZOTUs are labeled as <prefix><idx> where idx is determined by total counts across
    all samples, starting from <prefix>1 (e.g. 16S_ZOTU1).
    """
    if prefix is None:
        prefix = "ZOTU"
    # Load JSON data from file
    with open(input_json) as json_file:
        data = json.load(json_file)

    # Extract unique sequences and calculate their total counts, keeping track of "#UNKNOWN" separately
    sequence_counts = defaultdict(int)
    key_sequence_counts = defaultdict(lambda: defaultdict(int))
    unknown_counts = defaultdict(int)
    for key, value in data.items():
        for sequence, count in zip(value["zotus"], value["counts"]):
            if sequence == "#UNKNOWN":
                unknown_counts[key] = count
            else:
                sequence_counts[sequence] += count
                key_sequence_counts[key][sequence] = count

    # Sort sequences by total counts from high to low
    sorted_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)

    # Create a list of SeqRecord objects with names ZOTU1, ZOTU2, ..., and total counts as description
    seq_records = []
    sequence_to_zotu = {}
    for i, (sequence, count) in enumerate(sorted_sequences, start=1):
        zotu_name = f"{prefix}{i}"
        seq_records.append(
            SeqRecord(Seq(sequence), id=zotu_name, description=f"total_counts={count}")
        )
        sequence_to_zotu[sequence] = zotu_name

    # Write to FASTA file
    with open(output_fasta, "w") as fasta_file:
        SeqIO.write(seq_records, fasta_file, "fasta")

    # Create a DataFrame with ZOTU names as index and counts as values
    df_data = defaultdict(list)
    index = []

    for sequence, zotu_name in sequence_to_zotu.items():
        index.append(zotu_name)
        for key in key_sequence_counts:
            df_data[key].append(key_sequence_counts[key].get(sequence, 0))

    # Add #UNKNOWN to the DataFrame as the last row
    # if "#UNKNOWN" in unknown_count:
    if unknown_counts:
        index.append(f"{prefix}_UNKNOWN")
        for key in key_sequence_counts:
            df_data[key].append(unknown_counts[key])
        # df_data[f"{prefix}_UNKNOWN"] = df_data.pop("#UNKNOWN")

    df = pd.DataFrame(df_data, index=index)
    df.index.name = "#OTU_ID"

    # Save the DataFrame to a CSV file (optional)
    write_table(df, output_count)


def cluster_unoise3(
    uniq_fasta: str, out_fasta: str, minsize: int, relabel: str = "ZOTU"
):
    if os.path.isdir(out_fasta):
        raise ValueError(f"Output fasta file {out_fasta} is a directory")
    output_dir = os.path.dirname(out_fasta)
    os.makedirs(output_dir, exist_ok=True)
    # zotu_init = os.path.join(output_dir, f"_{out_fasta}")
    # vsearch does not have `--otutab`, `--otutab_stats`, `--calc_distmx`, and
    # `--cluster_aggd`. As suggested on the github issue
    #  (https://github.com/torognes/vsearch/issues/392), we use `--search_global`
    # to replace `--otutab`, but the other commands are not implemented.
    # Besides, in usearch chimera is removed as part of the clustering step, but
    # in vsearch chimera removal is a separate step.
    unoise3_proc = subprocess.Popen(
        [
            "vsearch",
            "--cluster_unoise",
            uniq_fasta,
            "--minsize",
            str(minsize),
            "--centroids",
            # f"{output_dir}/temp.fa",
            "-",
            # "--strand",
            # "both",
        ],
        stdout=subprocess.PIPE,
    )
    uchime3_proc = subprocess.Popen(
        [
            "vsearch",
            "--uchime3_denovo",
            "-",
            "--nonchimeras",
            out_fasta,
            "--relabel",
            relabel,
        ],
        stdin=unoise3_proc.stdout,
    )
    uchime3_proc.wait()
    # os.remove(zotu_init)


def _fq2fa(input_fastq: str, output_fasta: str) -> None:
    """Convert fastq (could be gzipped) to output fasta (gzipped) using seqtk"""
    try:
        # Constructing the seqtk command
        cmd = f"seqtk seq -a {input_fastq} | pigz -f > {output_fasta}"

        # Executing the command
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Error executing command: {e.stderr}")


def _cluster_unoise3(
    input_fastq: str,
    db_fasta: str,
    minsize: int,
    id_: float,
    output_dir: str | None,
    num_threads: int,
    backend: str = "vsearch",
) -> None:
    if output_dir is None:
        output_dir = os.path.dirname(input_fastq)
    os.makedirs(output_dir, exist_ok=True)
    if backend == "usearch":
        subprocess.run(
            [
                "usearch11",
                "-unoise3",
                db_fasta,
                "-zotus",
                f"{output_dir}/unoise3_zotu.fa",
                "-tabbedout",
                f"{output_dir}/unoise3_zotu.tsv",
            ],
            check=True,
        )
        subprocess.run(
            [
                "usearch11",
                "-fastx_relabel",
                f"{output_dir}/unoise3_zotu.fa",
                "-prefix",
                "ZOTU",
                "-fastaout",
                f"{output_dir}/unoise3_zotu_relabel.fa",
                "-keep_annots",
            ],
            check=True,
        )
        os.rename(
            f"{output_dir}/unoise3_zotu_relabel.fa", f"{output_dir}/unoise3_zotu.fa"
        )
        subprocess.run(
            [
                "usearch10",
                "-otutab",
                input_fastq,
                "-db",
                f"{output_dir}/unoise3_zotu.fa",
                "-otutabout",
                f"{output_dir}/unoise3_zotu.tsv",
                "-biomout",
                f"{output_dir}/unoise3_zotu.biom",
                "-mapout",
                f"{output_dir}/unoise3_zotu.map",
                "-notmatched",
                f"{output_dir}/unoise3_zotu_notmatched.fa",
                "-dbmatched",
                f"{output_dir}/unoise3_zotu_dbmatched.fa",
                "-sizeout",
                "-strand",
                "both",
                "-threads",
                str(num_threads),
                "-sample_delim",
                ".",
            ],
            check=True,
        )
        subprocess.run(
            [
                "usearch10",
                "-otutab_stats",
                f"{output_dir}/unoise3_zotu.tsv",
                "-output",
                f"{output_dir}/unoise3_zotu_report.txt",
            ],
            check=True,
        )
        subprocess.run(
            [
                "usearch11",
                "-calc_distmx",
                f"{output_dir}/unoise3_zotu.fa",
                "-tabbedout",
                f"{output_dir}/unoise3_zotu.dist",
                "-maxdist",
                "0.2",
                "-termdist",
                "0.3",
            ],
            check=True,
        )
        subprocess.run(
            [
                "usearch11",
                "-cluster_aggd",
                f"{output_dir}/unoise3_zotu.dist",
                "-treeout",
                f"{output_dir}/unoise3_zotu.tree",
                "-clusterout",
                f"{output_dir}/unoise3_zotu.clust",
                "-id",
                "0.8",
                "-linkage",
                "min",
            ],
            check=True,
        )
    elif backend == "vsearch":
        # vsearch does not have `--otutab`, `--otutab_stats`, `--calc_distmx`, and
        # `--cluster_aggd`. As suggested on the github issue
        #  (https://github.com/torognes/vsearch/issues/392), we use `--search_global`
        # to replace `--otutab`, but the other commands are not implemented.
        # Besides, in usearch chimera is removed as part of the clustering step, but
        # in vsearch chimera removal is a separate step.
        subprocess.run(
            [
                "vsearch",
                "--cluster_unoise",
                db_fasta,
                "--minsize",
                str(minsize),
                "--centroids",
                f"{output_dir}/temp.fa",
                "--strand",
                "both",
            ],
            check=True,
        )
        subprocess.run(
            [
                "vsearch",
                "--uchime3_denovo",
                f"{output_dir}/temp.fa",
                "--nonchimeras",
                f"{output_dir}/unoise3_zotu.fa",
                "--relabel",
                "ZOTU",
            ],
        )
        os.remove(f"{output_dir}/temp.fa")
        # convert input fastq to fasta
        with open(f"{os.path.splitext(input_fastq)[0]}.fa", "w") as f:
            subprocess.run(["seqtk", "seq", "-A", input_fastq], stdout=f, check=True)
        subprocess.run(
            [
                "vsearch",
                "--usearch_global",
                f"{os.path.splitext(input_fastq)[0]}.fa",
                "--db",
                f"{output_dir}/unoise3_zotu.fa",
                "--strand",
                "plus",
                "--id",
                str(id_),
                "--otutabout",
                f"{output_dir}/unoise3_zotu.tsv",
                "--biomout",
                f"{output_dir}/unoise3_zotu.biom",
                "--alnout",
                f"{output_dir}/unoise3_zotu.aln",
                "--dbmatched",
                f"{output_dir}/unoise3_zotu_dbmatched.fa",
                "--dbnotmatched",
                f"{output_dir}/unoise3_zotu_dbnotmatched.fa",
                "--sizeout",
                "--threads",
                str(num_threads),
            ]
        )


def tax_nbc(
    input_fasta: str,
    db_fasta: str | None,
    output_path: str,
    num_threads: str,
    backend: str = "rdp_classifier",
) -> None:
    if backend == "usearch":
        if db_fasta is None:
            raise ValueError("db_fasta is required when using usearch")
        _ = subprocess.run(
            [
                "usearch11",
                "-nbc_tax",
                input_fasta,
                "-db",
                db_fasta,
                "-strand",
                "both",
                "-tabbedout",
                output_path,
                "-threads",
                str(num_threads),
            ],
            check=True,
        )
    elif backend == "rdp_classifier":
        if db_fasta is not None:
            print("WARNING: db_fasta is ignored when using rdp_classifier")
        subprocess.run(["classifier", "-o", output_path, input_fasta], check=True)
    else:
        raise ValueError(f"Invalid backend: {backend}")


def tax_sintax(
    input_fastq: str,
    db_udb: str,
    output_path: str,
    num_threads: str,
    backend: str = "vsearch",
) -> None:
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    if backend == "usearch":
        subprocess.run(
            [
                "usearch10",
                "-sintax",
                input_fastq,
                "-db",
                db_udb,
                "-strand",
                "both",
                "-sintax_cutoff",
                "0.8",
                "-threads",
                str(num_threads),
                "-tabbedout",
                output_path,
            ],
            check=True,
        )
    elif backend == "vsearch":
        subprocess.run(
            [
                "vsearch",
                "--sintax",
                input_fastq,
                "--db",
                db_udb,
                "--strand",
                "both",
                "--sintax_cutoff",
                "0.8",
                "--threads",
                str(num_threads),
                "--tabbedout",
                output_path,
            ],
            check=True,
        )


def main():
    parser = argparse.ArgumentParser(description="Perform QC on FASTQ files")
    subparsers = parser.add_subparsers(dest="subcommand")

    merge_pairs_parser = subparsers.add_parser(
        "merge_pairs", help="Merge paired-end reads"
    )
    merge_pairs_parser.add_argument(
        "-i", "--fastq_dir", help="Input directory containing paired-end fastq files"
    )
    merge_pairs_parser.add_argument(
        "-o", "--fastq_out", help="Output path for merged fastq file"
    )
    merge_pairs_parser.add_argument(
        "-t", "--num_threads", type=int, default=16, help="Number of threads to use"
    )

    db_construct_parser = subparsers.add_parser(
        "db_construct", help="Construct a database from a FASTA file"
    )
    db_construct_parser.add_argument(
        "-i", "--input_fastq", help="Input FASTA file to construct the database from"
    )
    db_construct_parser.add_argument(
        "-o",
        "--output_fasta",
        help="Output path for database FASTA file",
        type=str,
        default=None,
    )
    db_construct_parser.add_argument(
        "-t", "--num_threads", type=int, default=16, help="Number of threads to use"
    )

    cluster_uparse_parser = subparsers.add_parser(
        "cluster_uparse", help="Cluster reads using UPARSE"
    )
    cluster_uparse_parser.add_argument(
        "-i", "--input_fastq", help="Input FASTQ file to cluster", type=str
    )
    cluster_uparse_parser.add_argument(
        "-d", "--db_fasta", help="Input FASTA file to cluster against", type=str
    )
    cluster_uparse_parser.add_argument(
        "-o",
        "--output_dir",
        help="Output directory for UPARSE results",
        type=str,
        default=None,
    )
    cluster_uparse_parser.add_argument(
        "-t", "--num_threads", type=int, default=16, help="Number of threads to use"
    )

    cluster_unoise3_parser = subparsers.add_parser(
        "cluster_unoise3", help="Cluster reads using UNOISE3"
    )
    cluster_unoise3_parser.add_argument(
        "-i", "--uniq_fasta", help="Input FASTQ file to cluster", type=str
    )
    cluster_unoise3_parser.add_argument(
        "-m", "--minsize", type=int, default=8, help="Minimum cluster size"
    )
    cluster_unoise3_parser.add_argument(
        "-o", "--out_fasta", help="Output path for ZOTU fasta file"
    )
    unoise3_parser = subparsers.add_parser(
        "unoise3", help="Cluster reads using UNOISE3"
    )
    unoise3_parser.add_argument(
        "-i",
        "--input_fastq",
        help="Input FASTQ file to cluster",
        type=str,
        required=True,
    )
    unoise3_parser.add_argument(
        "-o",
        "--output_fasta",
        help="Output path for ZOTU fasta file",
        type=str,
        required=True,
    )
    unoise3_parser.add_argument(
        "-m", "--minsize", type=int, default=8, help="Minimum cluster size"
    )
    unoise3_parser.add_argument("-a", "--alpha", type=float, default=2.0)

    unoise3_parser.add_argument(
        "--min_length",
        type=int,
        default=100,
        help="Minimum length of sequences to cluster",
    )
    unoise3_parser.add_argument(
        "-l", "--relabel_prefix", type=str, default=None, help="Prefix for ZOTU labels"
    )
    unoise3_parser.add_argument(
        "-t", "--num_threads", type=int, default=4, help="Number of threads to use"
    )
    search_global_parser = subparsers.add_parser(
        "search_global", help="Search reads against a database"
    )
    search_global_parser.add_argument(
        "-i", "--input_fastq", help="Input FASTQ file to search", type=str
    )
    search_global_parser.add_argument(
        "-d", "--zotu_fasta", help="Input FASTA file to search against", type=str
    )
    search_global_parser.add_argument(
        "-o",
        "--output_tsv",
        help="Output path for ZOTU tsv file",
        type=str,
        default=None,
    )
    search_global_parser.add_argument(
        "--id", type=float, default=0.97, help="Minimum cluster identity"
    )
    search_global_parser.add_argument(
        "--unknown_name",
        type=str,
        default="ZOTU_UNKNOWN",
        help="Name for unknown ZOTUs",
    )
    search_global_parser.add_argument(
        "-t", "--num_threads", type=int, default=8, help="Number of threads to use"
    )
    workflow_per_sample_parser = subparsers.add_parser(
        "workflow_per_sample", help="Run UNOISE3 workflow per sample"
    )
    workflow_per_sample_parser.add_argument(
        "-i", "--input_fastq", help="Input FASTQ file to cluster", type=str
    )
    workflow_per_sample_parser.add_argument(
        "-o", "--output_path", help="Output path", type=str
    )
    workflow_per_sample_parser.add_argument(
        "-m", "--min_size", type=int, default=8, help="Minimum cluster size"
    )
    workflow_per_sample_parser.add_argument("-a", "--alpha", type=float, default=2.0)
    workflow_per_sample_parser.add_argument(
        "-l",
        "--relabel_prefix",
        type=str,
        default="ZOTU",
        help="Prefix for ZOTU labels",
    )
    workflow_per_sample_parser.add_argument(
        "-t", "--num_threads", type=int, default=16, help="Number of threads to use"
    )
    workflow_per_sample_parser.add_argument(
        "--search",
        action="store_true",
        help="Search reads against ZOTU database",
    )
    aggregate_samples_parser = subparsers.add_parser(
        "aggregate_samples", help="Aggregate ZOTU counts across samples"
    )
    aggregate_samples_parser.add_argument(
        "-i", "--input_json", help="Input JSON file with ZOTU counts", type=str
    )
    aggregate_samples_parser.add_argument(
        "-of", "--output_fasta", help="Output path for aggregated ZOTU fasta", type=str
    )
    aggregate_samples_parser.add_argument(
        "-oc", "--output_count", help="Output path for aggregated ZOTU counts", type=str
    )
    aggregate_samples_parser.add_argument(
        "-p",
        "--prefix",
        help="Prefix for ZOTU labels",
        type=str,
        default=None,
    )

    tax_nbc_parser = subparsers.add_parser(
        "tax_nbc", help="Taxonomy classification using NBC"
    )
    tax_nbc_parser.add_argument(
        "-i", "--input_fasta", help="Input FASTA file to classify", type=str
    )
    tax_nbc_parser.add_argument(
        "-d", "--db_fasta", help="Input FASTA file to classify against", type=str
    )
    tax_nbc_parser.add_argument(
        "-o", "--output_path", help="Output path for NBC results", type=str
    )
    tax_nbc_parser.add_argument(
        "-t", "--num_threads", type=int, default=16, help="Number of threads to use"
    )

    tax_sintax_parser = subparsers.add_parser(
        "tax_sintax", help="Taxonomy classification using SINTAX"
    )
    tax_sintax_parser.add_argument(
        "-i", "--input_fastq", help="Input FASTQ file to classify", type=str
    )
    tax_sintax_parser.add_argument(
        "-d", "--db_udb", help="Input udb file to classify against", type=str
    )
    tax_sintax_parser.add_argument(
        "-o", "--output_path", help="Output path for SINTAX results", type=str
    )
    tax_sintax_parser.add_argument(
        "-t", "--num_threads", type=int, default=16, help="Number of threads to use"
    )

    args = parser.parse_args()

    if args.subcommand == "db_construct":
        fastx_post_qc = qc(args.input_fastq, None, args.num_threads)
        db_construct(fastx_post_qc, args.output_fasta, args.num_threads)
    elif args.subcommand == "cluster_uparse":
        cluster_uparse(
            args.input_fastq, args.db_fasta, args.output_dir, args.num_threads
        )
    elif args.subcommand == "cluster_unoise3":
        cluster_unoise3(args.uniq_fasta, args.minsize, args.out_fasta)
    elif args.subcommand == "unoise3":
        unoise3(
            args.input_fastq,
            args.output_fasta,
            min_size=args.minsize,
            min_length=args.min_length,
            relabel_prefix=args.relabel_prefix,
            num_threads=args.num_threads,
        )
    elif args.subcommand == "search_global":
        search_global_add_unknown(
            args.input_fastq,
            args.zotu_fasta,
            args.output_tsv,
            id_=args.id,
            unknown_name=args.unknown_name,
            num_threads=args.num_threads,
        )
    elif args.subcommand == "workflow_per_sample":
        workflow_per_sample(
            args.input_fastq,
            args.output_path,
            min_size=args.min_size,
            alpha=args.alpha,
            prefix=args.relabel_prefix,
            num_threads=args.num_threads,
            search=args.search,
        )
    elif args.subcommand == "aggregate_samples":
        aggregate_samples(
            args.input_json, args.output_fasta, args.output_count, prefix=args.prefix
        )
    elif args.subcommand == "tax_nbc":
        tax_nbc(args.input_fasta, args.db_fasta, args.output_path, args.num_threads)
    elif args.subcommand == "tax_sintax":
        tax_sintax(args.input_fastq, args.db_udb, args.output_path, args.num_threads)
    else:
        raise ValueError(f"Invalid subcommand: {args.subcommand}")


if __name__ == "__main__":
    main()
