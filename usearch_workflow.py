#!python3

import argparse
import glob
import io
import os
import subprocess
import shutil

import numpy as np
import pandas as pd
from tqdm.auto import tqdm


def merge_pairs(
    fastq_dir: str, output_dir: str, num_threads: int, backend: str = "vsearch"
) -> None:
    # merge paired end reads
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
        for f in os.listdir(fastq_dir):
            if f.endswith(".fastq"):
                sample_name = f.split("_")[0]
                # use Popen for non-blocking execution
                outfile = open(f"{output_dir}/temp/{f}", "w")
                process = subprocess.Popen(
                    [
                        "seqtk",
                        "rename",
                        f"{fastq_dir}/{f}",
                        f"sample={sample_name} ",
                    ],
                    stdout=outfile,
                )
                fs.append(outfile)
                ps.append(process)
        [p.wait() for p in ps]
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
                f"{output_dir}/merged.fastq",
                "--fastaout",
                f"{output_dir}/merged.fa",
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
    input_fastq: str, num_splits: int, num_threads: int = 16, backend="vsearch"
) -> None:
    output_fastq = os.path.splitext(input_fastq)[0] + ".filtered.fa"

    if backend == "usearch":
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
            with open(output_fastq, "w") as outfile:
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
                    output_fastq,
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
                "--fastaout",
                output_fastq,
                "--fastq_maxee",
                "0.5",
                "--fastq_minlen",
                "200",
                "--fastq_maxns",
                "0",
                "--relabel",
                "filtered",
                "--threads",
                str(num_threads),
            ]
        )


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
    command = "seqkit stats input_fastqs/*_R1.fastq -T"
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
        raise ValueError(
            f"These samples in the fastq directory are not present in the metadata: {extra_samples}"
        )
    extra_samples = df.index.difference(df_rc.index)
    if len(extra_samples) > 0:
        print(
            f"WARNING: These samples in the metadata are not present in the fastq directory: {extra_samples}"
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
        read_count = min(row["read_count_subsample"], row["read_count"])
        num_subs = (
            num_subsamples if row["read_count_subsample"] < row["read_count"] else 1
        )
        del row["read_count_subsample"], row["read_count"]
        ps = []
        for iteration in range(1, num_subs + 1):
            # Increment the seed value
            seed += 1
            # Subsample R1 and R1 and output to the subsampled directory
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
            subsample_metadata.append(
                {
                    "sample": f"{sample_name}-subsampled-{iteration}",
                    "original_sample": sample_name,
                    "read_count": read_count,
                }
                | row.to_dict()
            )
        [p.wait() for p in ps]
    pd.DataFrame(subsample_metadata).to_csv(output_metadata_path, sep="\t", index=False)


def db_construct(input_fasta, num_threads: int, backend: str = "vsearch"):
    output_dir, base_name = os.path.split(input_fasta)
    base_name_noext = os.path.splitext(base_name)[0]

    if backend == "usearch":
        subprocess.run(
            [
                "usearch11",
                "-fastx_uniques",
                input_fasta,
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
                input_fasta,
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
                "--fastx_unique",
                input_fasta,
                "--fastaout",
                f"{output_dir}/{base_name_noext}.uniq.fa",
                "--sizeout",
            ]
        )
        subprocess.run(
            [
                "vsearch",
                "--sortbysize",
                f"{output_dir}/{base_name_noext}.uniq.fa",
                "--output",
                f"{output_dir}/{base_name_noext}.sorted.fa",
                # subsequent clustering/denoising will apply threshold, so not necessary to do it here.
                # "--minsize",
                # "2",
            ]
        )
        subprocess.run(
            [
                "vsearch",
                "--search_exact",
                input_fasta,
                "--db",
                f"{output_dir}/{base_name_noext}.sorted.fa",
                "--strand",
                "both",
                "--notmatched",
                f"{output_dir}/{base_name_noext}.unmatched.fa",
                "--threads",
                str(num_threads),
            ]
        )


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


def cluster_unoise3(
    input_fastq: str,
    db_fasta: str,
    output_dir: str,
    num_threads: int,
    backend="vsearch",
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
        # instead, but the other commands are not implemented.
        # Besides, in usearch chimera is removed as part of the clustering step, but
        # in vsearch chimera removal is a separate step.
        subprocess.run(
            [
                "vsearch",
                "--cluster_unoise",
                db_fasta,
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
                "0.97",
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
    db_fasta: str,
    output_path: str,
    num_threads: str,
    backend: str = "rdp_classifier",
) -> None:
    if backend == "usearch":
        subprocess.run(
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

    # add depth subcommand
    add_depth_to_metadata_parser = subparsers.add_parser(
        "add_depth_to_metadata", help="Add read depth to metadata file"
    )
    add_depth_to_metadata_parser.add_argument(
        "-i", "--fastq_dir", help="Input directory containing paired-end fastq files"
    )
    add_depth_to_metadata_parser.add_argument(
        "-m", "--metadata_path", help="Path to metadata file"
    )

    # Subsample subcommand
    subsample_parser = subparsers.add_parser(
        "subsample", help="Subsample input FASTQ files"
    )
    subsample_parser.add_argument(
        "-i", "--fastq_dir", help="Input directory containing paired-end fastq files"
    )
    subsample_parser.add_argument("-m", "--metadata_path", help="Path to metadata file")
    subsample_parser.add_argument(
        "-o", "--output_dir", help="Output directory for subsampled fastq files"
    )
    subsample_parser.add_argument(
        "-om", "--output_metadata_path", help="Output metadata file"
    )
    subsample_parser.add_argument(
        "-n", "--num_subsamples", type=int, help="Number of subsamples to generate"
    )
    subsample_parser.add_argument(
        "--mode",
        choices=["reference", "min"],
        default="reference",
        help="Subsampling mode",
    )
    subsample_parser.add_argument(
        "-s", "--seed", type=int, default=100, help="Random seed"
    )

    merge_pairs_parser = subparsers.add_parser(
        "merge_pairs", help="Merge paired-end reads"
    )
    merge_pairs_parser.add_argument(
        "-i", "--fastq_dir", help="Input directory containing paired-end fastq files"
    )
    merge_pairs_parser.add_argument(
        "-o", "--output_dir", help="Output directory for merged fastq file"
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
        "-s",
        "--num_splits",
        type=int,
        help="Number of splits to split the FASTQ file into",
        default=1,
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
        "-i", "--input_fastq", help="Input FASTQ file to cluster", type=str
    )
    cluster_unoise3_parser.add_argument(
        "-d", "--db_fasta", help="Input FASTA file to cluster against", type=str
    )
    cluster_unoise3_parser.add_argument(
        "-o",
        "--output_dir",
        help="Output directory for UNOISE3 results",
        type=str,
        default=None,
    )
    cluster_unoise3_parser.add_argument(
        "-t", "--num_threads", type=int, default=16, help="Number of threads to use"
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

    if args.subcommand == "add_depth_to_metadata":
        add_depth_to_metadata(args.fastq_dir, args.metadata_path)
    elif args.subcommand == "subsample":
        subsample(
            args.fastq_dir,
            args.metadata_path,
            args.output_dir,
            args.output_metadata_path,
            args.num_subsamples,
            args.mode,
            args.seed,
        )
    elif args.subcommand == "merge_pairs":
        merge_pairs(args.fastq_dir, args.output_dir, args.num_threads)
    # elif args.subcommand == "qc":
    #     qc(args.input_fastq, args.num_splits, args.num_threads)
    elif args.subcommand == "db_construct":
        qc(args.input_fastq, args.num_splits, args.num_threads)
        db_construct(
            os.path.splitext(args.input_fastq)[0] + ".filtered.fa", args.num_threads
        )
    elif args.subcommand == "cluster_uparse":
        cluster_uparse(
            args.input_fastq, args.db_fasta, args.output_dir, args.num_threads
        )
    elif args.subcommand == "cluster_unoise3":
        cluster_unoise3(
            args.input_fastq, args.db_fasta, args.output_dir, args.num_threads
        )
    elif args.subcommand == "tax_nbc":
        tax_nbc(args.input_fasta, args.db_fasta, args.output_path, args.num_threads)
    elif args.subcommand == "tax_sintax":
        tax_sintax(args.input_fastq, args.db_udb, args.output_path, args.num_threads)


if __name__ == "__main__":
    main()
