#!/usr/bin/env python

import glob
import os
import shutil
import subprocess
import argparse
from pathlib import Path
import resource

from Bio.Seq import Seq

from masato.utils import cat_fastq, cat_fastq_se, print_command


def get_rc(seq: str) -> str:
    """Return the reverse complement of the given sequence."""
    return str(Seq(seq).reverse_complement())


TRUSEQ_READ1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
TRUSEQ_READ2 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
PRIMER_ITS_5 = "GTACCTGCGGARGGATCA"
PRIMER_ITS_7 = "ATGAGATCCRTTGYTRAAAGTT"
PRIMER_16S_5 = "GTGTGYCAGCMGCCGCGGTAA"
PRIMER_16S_7 = "GGACTACNVGGGTWTCTAAT"
PRIMER_16S_V3V4_341F = "CCTACGGGNGGCWGCAG"
PRIMER_16S_V3V4_805R = "GACTACHVGGGTATCTAATCC"
PRIMER_TN5 = "AGATGTGTATAAGAGACAG"
# PRIMER_MAPS_0 = "GCGAGACACTGCTGAGCG"
# PRIMER_MAPS_1 = "GTCACTCCAGCCTCGTCG"
# PRIMER_MAPS_2 = "AGGCAGCGTCTGTAGCGA"
PRIMER_MAPS_1 = "GCGAGACACTGCTGAGCG"
PRIMER_MAPS_2 = "GTCACTCCAGCCTCGTCG"
PRIMER_MAPS_3 = "GTGTTCGTCGGCAGCGTC" + PRIMER_TN5
PRIMER_MAPS_r = "GTCTCGTGGGCTCGG" + PRIMER_TN5
ECREC_SPACER0 = "GAGCACAAATATCATCGCTCAAACC"
ECREC_DR = "GTGTTCCCCGCGCCAGCGGGGATAAACC"
ECREC_LEADER = "CTGGCTTAAAAAATCATTAATTAATAATAGGTTATGTTTAGA"
RECORDING_PRIMER_3 = "AGATCGGAAGAGCACACGTCTGA"
RECORDING_PRIMER_5 = "CCTACACGACGCTCTTCCGATCT"


def set_nofile_limit():
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    resource.setrlimit(resource.RLIMIT_NOFILE, (100000, min(100000, hard)))


def run_trim_galore(input_dir, output_dir, pair):
    """Run trim_galore on the given pair of files."""
    r1_file, r2_file = pair
    cmd = [
        "trim_galore",
        "--paired",
        os.path.join(input_dir, r1_file),
        os.path.join(input_dir, r2_file),
        "-a",
        "GGCTGACTGACT",
        "-o",
        output_dir,
        "-j",
        "4",
    ]
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def run_trimmomatic(input_fastq_1: str, input_fastq_2: str, output_dir: str) -> None:
    """Run trimmomatic on the given pair of files using subprocess, following this example:
    trimmomatic-0.30.jar PE s_1_1_sequence.txt.gz s_1_2_sequence.txt.gz
    lane1_forward_paired.fq.gz lane1_forward_unpaired.fq.gz lane1_reverse_paired.fq.gz
    lane1_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3
    TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
    output_fastq_1 = os.path.join(output_dir, os.path.basename(input_fastq_1))
    output_fastq_2 = os.path.join(output_dir, os.path.basename(input_fastq_2))
    output_fastq_1_unpaired = output_fastq_1.replace(
        ".fastq", ".trimmomatic.unpaired.fastq"
    )
    output_fastq_2_unpaired = output_fastq_2.replace(
        ".fastq", ".trimmomatic.unpaired.fastq"
    )
    output_trim_log = os.path.join(output_dir, "trimmomatic.log")
    output_summary = os.path.join(output_dir, "trimmomatic.summary")
    os.makedirs(output_dir, exist_ok=True)
    subprocess.run(
        [
            "trimmomatic",
            "PE",
            "-threads",
            "8",
            "-trimlog",
            output_trim_log,
            "-summary",
            output_summary,
            input_fastq_1,
            input_fastq_2,
            output_fastq_1,
            output_fastq_1_unpaired,
            output_fastq_2,
            output_fastq_2_unpaired,
            # "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10",
            "LEADING:3",
            "TRAILING:3",
            "SLIDINGWINDOW:4:15",
            "MINLEN:36",
        ]
    )


def run_fastp(
    input_fastq_1: str, input_fastq_2: str, output_dir: str, stdout: bool = False
) -> subprocess.Popen[bytes]:
    output_fastq_1 = os.path.join(output_dir, os.path.basename(input_fastq_1))
    output_fastq_2 = os.path.join(output_dir, os.path.basename(input_fastq_2))
    output_fastq_1_unpaired = output_fastq_1.replace(".fastq", ".fastp.unpaired.fastq")
    output_fastq_2_unpaired = output_fastq_2.replace(".fastq", ".fastp.unpaired.fastq")
    failed_out = output_fastq_2.replace(".fastq", ".fastp.failed.fastq")
    report_out = os.path.join(output_dir, "fastp.report.json")
    os.makedirs(output_dir, exist_ok=True)
    fastp_proc = subprocess.Popen(
        [
            "fastp",
            "-i",
            input_fastq_1,
            "-I",
            input_fastq_2,
            "--unpaired1",
            output_fastq_1_unpaired,
            "--unpaired2",
            output_fastq_2_unpaired,
            "--failed_out",
            failed_out,
            # cut 5'
            "-5",
            "--cut_front_window_size",
            "1",
            "--cut_front_mean_quality",
            "3",
            # cut 3'
            "-3",
            "--cut_tail_window_size",
            "1",
            "--cut_tail_mean_quality",
            "3",
            # cut right
            "-r",
            "--cut_right_window_size",
            "4",
            "--cut_right_mean_quality",
            "15",
            "--length_required",
            "36",
            "--json",
            report_out,
            "-w",
            "8",
        ]
        + (
            [
                "-o",
                output_fastq_1,
                "-O",
                output_fastq_2,
            ]
            if stdout is False
            else ["--stdout"]
        ),
        stdout=None if stdout is False else subprocess.PIPE,
    )
    return fastp_proc


def run_fastp_cutadapt(
    input_fastq_1: str,
    input_fastq_2: str,
    output_dir: str,
    amplicon_type: str = "ITS",
):
    """Quality trimming and filtering by fastp, then adapter trimming by cutadapt.
    For read 1, trim RC of i7 adapter from the 3' end; For read 2, trim RC of i5 adapter
        also from the 3' end. Just trim and save report, no need to separate trimmed
        and untrimmed reads or discard reads.
    """
    output_fastq_1: str = os.path.join(output_dir, os.path.basename(input_fastq_1))
    output_fastq_2: str = os.path.join(output_dir, os.path.basename(input_fastq_2))
    report_report: str = os.path.join(output_dir, "cutadapt.report.json")
    if amplicon_type not in ["ITS", "16S"]:
        raise ValueError("amplicon_type must be either ITS or 16S.")

    fastp_proc = run_fastp(input_fastq_1, input_fastq_2, output_dir, stdout=True)
    _ = subprocess.run(
        [
            "cutadapt",
            "-a",
            get_rc(PRIMER_ITS_7 if amplicon_type == "ITS" else PRIMER_16S_7),
            "-A",
            get_rc(PRIMER_ITS_5 if amplicon_type == "ITS" else PRIMER_16S_5),
            "-o",
            output_fastq_1,
            "-p",
            output_fastq_2,
            "--json",
            report_report,
            "--interleaved",
            "-",
        ],
        stdin=fastp_proc.stdout,
    )
    fastp_proc.wait()


def simple_preprocess(fastq_dir: str, output_fastq: str) -> None:
    """Given a directory of fastq files, perform renaming, quality trimming, adapter
    trimming and length filtering using just python and fastp. fastp takes interleaved
    input from stdin, to which we write the renamed and interleaved paired-end reads.

    Note:
    - Adapter sequences are optional since fastp can trim paired-end reads by overlap
        analysis, which is more versatile and robust.
    - fastp can also merge pairs, but we don't do it here but with the following ZOTU
        pipeline, just to be more coherent.
    """
    output_dir = os.path.dirname(output_fastq)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "fastp"), exist_ok=True)
    fastp_proc = subprocess.Popen(
        [
            "fastp",
            # "-i",
            # os.path.join(fastq_dir, "*_R1_001.fastq.gz"),
            # "-I",
            # os.path.join(fastq_dir, "*_R2_001.fastq.gz"),
            "--stdin",
            "--interleaved_in",
            output_fastq,
            "-w",  # number of threads
            "8",
            # don't trim, so that amplicons have uniform length, as suggested by Edgar
            "--length_required",
            "200",
            # if you don't want to merge, modify the following arguments
            # for the meaning of each of the following arguments, see fastp github README
            "--merge",
            "--merged_out",
            output_fastq,
            "--unpaired1",
            os.path.join(output_dir, "fastp", "unpaired_R1.fastq.gz"),
            "--unpaired2",
            os.path.join(output_dir, "fastp", "unpaired_R2.fastq.gz"),
            "--failed_out",
            os.path.join(output_dir, "fastp", "failed.fastq.gz"),
            "--out1",
            os.path.join(output_dir, "fastp", "unmerged_R1.fastq.gz"),
            "--out2",
            os.path.join(output_dir, "fastp", "unmerged_R2.fastq.gz"),
            "--html",  # report
            os.path.join(output_dir, "fastp", "report.html"),
            "--json",
            os.path.join(output_dir, "fastp", "report.json"),
        ],
        stdin=subprocess.PIPE,
        # stdout=subprocess.PIPE,
    )
    cat_fastq(fastq_dir, output_fp_r1=fastp_proc.stdin, output_fp_r2=fastp_proc.stdin)
    fastp_proc.stdin.close()
    fastp_proc.wait()


def pseudo_merge(
    fastq_dir: str,
    output_fastq: str,
    threads: int,
    min_l1: int,
    first_k1: int,
    min_l2: int,
    first_k2: int,
    cutadapt_error: float = 0.1,
) -> None:
    """
    Process all FASTQ files in a directory by performing quality trimming, adapter trimming,
    and merging paired-end reads.

    Args:
        fastq_dir (str): Directory containing FASTQ files.
        output_dir (str): Directory to save processed FASTQ files.
        threads (int): Number of threads to use.
        min_l1 (int): Minimum length for R1 reads after trimming.
        first_k1 (int): Target length for R1 reads after trimming.
        min_l2 (int): Minimum length for R2 reads after trimming.
        first_k2 (int): Target length for R2 reads after trimming.
        cutadapt_error (float): Error tolerance for cutadapt trimming.
    """
    output_dir = os.path.dirname(output_fastq)
    os.makedirs(output_dir, exist_ok=True)
    log_dir = Path(output_dir) / "trim_logs"
    os.makedirs(log_dir, exist_ok=True)

    fastp_html_r1 = log_dir / "fastp.html"
    fastp_json_r1 = log_dir / "fastp.json"

    trimmed_r1 = Path(output_dir) / "trimmed_1.fastq.gz"
    trimmed_r2 = Path(output_dir) / "trimmed_2.fastq.gz"

    # Process R1 using cat_fastq_se() directly into fastp.
    # Unpaired or failed reads are discarded.
    fastp_args = [
        "fastp",
        "--stdin",
        "--interleaved_in",
        "--stdout",
        "--cut_right",
        "--html",
        str(fastp_html_r1),
        "--json",
        str(fastp_json_r1),
        "--thread",
        str(threads),
    ]
    print_command(fastp_args)
    fastp_proc = subprocess.Popen(
        fastp_args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )

    # too short reads are discarded.
    cutadapt_args = [
        "cutadapt",
        "--interleaved",
        "--minimum-length",
        f"{min_l1}:{min_l2}",
        "--length",
        f"{first_k1}",
        "-L",  # there is not a long version of this argument, but it is -l for R2. This is a 4.9 feature.
        f"{first_k2}",
        "--pair-filter",
        "both",
        "-e",
        str(cutadapt_error),
        "-o",
        str(trimmed_r1),
        "-p",
        "-",
        "--cores",
        str(threads),
        "-",
    ]
    print_command(cutadapt_args)
    cutadapt_proc = subprocess.Popen(
        cutadapt_args,
        stdin=fastp_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )

    # fastp_min_length_args = [
    #     "fastp",
    #     "--stdin",
    #     "--interleaved_in",
    #     "--out1",
    #     str(trimmed_r1),
    #     "--out2",  # out2 to stdout
    #     "/dev/stdout",
    #     "--max_len1",
    #     str(first_k1),
    #     "--max_len2",
    #     str(first_k2),
    #     "--thread",
    #     str(threads),
    # ]
    # print_command(fastp_min_length_args)
    # fastp_min_length_proc = subprocess.Popen(
    #     fastp_min_length_args,
    #     stdin=cutadapt_proc.stdout,
    #     stdout=subprocess.PIPE,
    #     stderr=subprocess.DEVNULL,
    # )

    seqkit_args = ["seqkit", "seq", "-rp", "-t", "DNA", "-o", str(trimmed_r2)]
    print_command(seqkit_args)
    seqkit_proc = subprocess.Popen(
        seqkit_args,
        stdin=cutadapt_proc.stdout,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    # Pipe R1 reads into fastp
    cat_fastq(fastq_dir, output_fp_r1=fastp_proc.stdin)

    fastp_proc.stdin.close()
    fastp_proc.stdout.close()
    cutadapt_proc.stdout.close()
    seqkit_proc.wait()
    # wait for all processes to finish
    # fastp_proc.wait()
    # cutadapt_proc.wait()
    # fastp_min_length_proc.wait()
    # seqkit_proc.wait()

    # Concatenate trimmed R1 and R2
    _merge_fastq_gz(str(trimmed_r1), str(trimmed_r2), output_fastq)
    os.remove(trimmed_r1)
    os.remove(trimmed_r2)


def rename_files_with_mmv(file_dir: str, patterns_file: str) -> None:
    # Step 1: Resolve absolute path of the patterns file
    patterns_file_abs = os.path.abspath(patterns_file)

    # Save the current working directory
    original_dir = os.getcwd()

    # Step 2: Change to the temporary directory
    os.chdir(file_dir)

    # Step 3: Execute the mmv command
    # with open(patterns_file_abs) as file:
    #     subprocess.run(["mmv", "-d"], stdin=file)
    with open(patterns_file_abs) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            old_name, new_name = line.split()
            subprocess.run(["mv", old_name, new_name])

    # Step 4: Change back to the original directory
    os.chdir(original_dir)


def get_min_overlap(adapter: str, frac: float = 0.8) -> int:
    return int(len(adapter) * frac)


def get_primer_set(name: str) -> tuple[str, str | None] | dict[str, str] | str:
    if name == "its":
        return (
            f"^{PRIMER_ITS_5};required...{get_rc(PRIMER_ITS_7)};optional",
            f"^{PRIMER_ITS_7};required...{get_rc(PRIMER_ITS_5)};optional",
        )
    elif name == "16s":
        return (
            f"^{PRIMER_16S_5};required...{get_rc(PRIMER_16S_7)};optional",
            f"^{PRIMER_16S_7};required...{get_rc(PRIMER_16S_5)};optional",
        )
    elif name == "16s_strange":
        primer_5 = f"{PRIMER_16S_5}"
        primer_7 = f"{PRIMER_16S_7}"
        return (
            f"{primer_5};required...{get_rc(primer_5)};optional",
            f"{primer_7};required...{get_rc(primer_7)};optional",
        )
    elif name == "its_strange":
        primer_5 = f"{PRIMER_ITS_5}"
        primer_7 = f"{PRIMER_ITS_7}"
        return (
            f"{primer_5};required...{get_rc(primer_5)};optional",
            f"{primer_7};required...{get_rc(primer_7)};optional",
        )
    elif name == "its_3":
        return get_rc(PRIMER_ITS_7), get_rc(PRIMER_ITS_5)
    elif name == "16s_3":
        return get_rc(PRIMER_16S_7), get_rc(PRIMER_16S_5)
    # elif name == "maps_0":
    #     # optional
    #     return f"{PRIMER_MAPS_0}...{get_rc(PRIMER_MAPS_r)}", None
    elif name == "maps_round1":
        # anchor
        return f"^{PRIMER_MAPS_1}", None
    elif name == "maps_round2":
        # anchor
        return f"^{PRIMER_MAPS_2}", None
    elif name == "maps_3":
        # anchor
        return f"^{PRIMER_MAPS_3}", None
    elif name == "maps_rand_hex_test":
        primer_maps = PRIMER_MAPS_1 + "N" * 6
        return (
            f"^{primer_maps};required...{get_rc(PRIMER_MAPS_r)};optional",
            get_rc(primer_maps),
        )
    elif name == "recording_adapter":
        return (
            f"{RECORDING_PRIMER_5};rightmost;min_overlap={get_min_overlap(RECORDING_PRIMER_5, frac=0.5)}"
            f"...{RECORDING_PRIMER_3};min_overlap={get_min_overlap(RECORDING_PRIMER_3, frac=0.5)}"
        )
    elif name == "recording_leader_dr":
        ret = ECREC_LEADER + ECREC_DR
        return ret + f";min_overlap={get_min_overlap(ret)}"
    elif name == "recording_spacer0":
        return (
            ECREC_SPACER0 + f";min_overlap={get_min_overlap(ECREC_SPACER0)};rightmost"
        )
    elif name == "recording_spacer0_rc":
        return get_rc(ECREC_SPACER0) + f";min_overlap={get_min_overlap(ECREC_SPACER0)}"
    elif name == "recording_dr":
        return ECREC_DR + f";min_overlap={get_min_overlap(ECREC_DR)}"
    elif name == "recording_dr_rc":
        return get_rc(ECREC_DR) + f"X;min_overlap={get_min_overlap(ECREC_DR)}"
    else:
        raise ValueError(f"Invalid primer set name: {name}")


def isolate_150_preprocess(
    fastq_dir: str,
    barcode_fwd_fasta: str,
    barcode_rev_fasta: str,
    rename_pattern: str,
    output_fastq: str,
    primer_set: str,
    first_k: int | None = None,
    min_length: int = 100,
    early_stop: bool = False,
) -> None:
    output_dir, output_f = os.path.split(output_fastq)
    output_dir_demux = os.path.join(output_dir, "demux")
    output_dir_demux_fail = os.path.join(output_dir, "demux_failed")
    output_dir_cutadapt = os.path.join(output_dir, "cutadapt")
    output_fastq_r1 = "_R1.".join(output_f.split(".", 1))
    output_fastq_r2 = "_R2.".join(output_f.split(".", 1))
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir_demux, exist_ok=True)
    os.makedirs(output_dir_demux_fail, exist_ok=True)
    os.makedirs(output_dir_cutadapt, exist_ok=True)

    # demultiplex with cutadapt
    proc_args = [
        "cutadapt",
        "-e",
        "0.15",
        "--no-indels",
        "-g",
        f"^file:{barcode_fwd_fasta}",
        "-G",
        f"^file:{barcode_rev_fasta}",
        "-o",
        os.path.join(output_dir_demux, "{name1}-{name2}_R1.fq.gz"),
        "-p",
        os.path.join(output_dir_demux, "{name1}-{name2}_R2.fq.gz"),
        "--cores",
        "1",
        "--interleaved",
        "-",
    ]
    print_command(proc_args)
    cutadapt_demux_proc = subprocess.Popen(
        proc_args,
        stdin=subprocess.PIPE,
        # silence cutadapt's output
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    cat_fastq(
        fastq_dir,
        output_fp_r1=cutadapt_demux_proc.stdin,
        output_fp_r2=cutadapt_demux_proc.stdin,
    )
    cutadapt_demux_proc.stdin.close()
    cutadapt_demux_proc.wait()
    rename_files_with_mmv(output_dir_demux, rename_pattern)

    # move all remaining files (whose name1 or name2 is "unknown") to output_dir/demux_failed
    demux_failed_files = list(
        set(
            glob.glob(os.path.join(output_dir_demux, "*-unknown_R1.fq.gz"))
            + glob.glob(os.path.join(output_dir_demux, "*-unknown_R2.fq.gz"))
            + glob.glob(os.path.join(output_dir_demux, "unknown-*_R1.fq.gz"))
            + glob.glob(os.path.join(output_dir_demux, "unknown-*_R2.fq.gz"))
        )
    )
    subprocess.run(
        [
            "mv",
            *demux_failed_files,
            os.path.join(output_dir, "demux_failed"),
        ]
    )

    # trim with cutadapt
    a, A = get_primer_set(primer_set)
    proc_args = (
        ["cutadapt", "-e", "0.15", "-a", a]
        + (["-A", A] if not early_stop else [])  # no trimming for read 2 if early stop
        + [
            "--minimum-length",
            str(min_length),
            "--pair-filter",
            "any",
            "-O",
            "16",
            "-o",
            os.path.join(output_dir_cutadapt, output_fastq_r1),
            "-p",
            os.path.join(output_dir_cutadapt, output_fastq_r2),
            "--untrimmed-output",
            os.path.join(output_dir_cutadapt, "untrimmed_1.fq.gz"),
            "--untrimmed-paired-output",
            os.path.join(output_dir_cutadapt, "untrimmed_2.fq.gz"),
            "--too-short-output",
            os.path.join(output_dir_cutadapt, "too_short_1.fq.gz"),
            "--too-short-paired-output",
            os.path.join(output_dir_cutadapt, "too_short_2.fq.gz"),
            "--cores",
            "4",
            "--interleaved",
            "-",
        ]
    )
    if first_k is not None:
        proc_args.extend(["-l", str(first_k)])
    print_command(proc_args)
    cutadapt_trim_proc = subprocess.Popen(
        proc_args,
        stdin=subprocess.PIPE,
        # silence cutadapt's output
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    cat_fastq(
        output_dir_demux,
        output_fp_r1=cutadapt_trim_proc.stdin,
        output_fp_r2=cutadapt_trim_proc.stdin,
        _have_sample_name=True,
    )
    cutadapt_trim_proc.stdin.close()
    cutadapt_trim_proc.wait()

    # copy merged_1 to output_fastq
    shutil.copy(os.path.join(output_dir_cutadapt, output_fastq_r1), output_fastq)
    shutil.rmtree(output_dir_demux)


def cutadapt_demux_se(fastq_path: str, output_dir: str, barcode_fastq: str):
    if not os.path.isfile(barcode_fastq):
        raise ValueError(f"{barcode_fastq} does not exist.")
    output_dir_cutadapt = os.path.join(output_dir, "cutadapt")
    output_dir_demux = os.path.join(output_dir, "demux")
    output_dir_demux_fail = os.path.join(output_dir, "demux_failed")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir_demux, exist_ok=True)
    os.makedirs(output_dir_demux_fail, exist_ok=True)
    os.makedirs(output_dir_cutadapt, exist_ok=True)
    # demultiplex with cutadapt
    proc_args = [
        "cutadapt",
        "-e",
        "0.2",
        "--no-indels",
        "-g",
        f"^file:{barcode_fastq}",
        "-o",
        os.path.join(output_dir_demux, "{name}.fq.gz"),
    ]
    proc_args.append(fastq_path)
    subprocess.run(proc_args)
    subprocess.run(
        [
            "mv",
            os.path.join(output_dir_demux, "unknown.fq.gz"),
            os.path.join(output_dir, "demux_failed"),
        ]
    )


def cutadapt_demux_merge_trim_se(
    fastq_path: str,
    output_fastq: str,
    primer_set: str,
    barcode_fastq: str,
    first_k: int | None = None,
    min_length: int | None = None,
    _r2: bool = False,
) -> None:
    """Demultiplex and merge single-end reads using cutadapt."""
    if not os.path.isfile(barcode_fastq):
        raise ValueError(f"{barcode_fastq} does not exist.")
    output_dir, output_f = os.path.split(output_fastq)
    output_dir_cutadapt = os.path.join(output_dir, "cutadapt")
    output_dir_demux = os.path.join(output_dir, "demux")
    output_dir_demux_fail = os.path.join(output_dir, "demux_failed")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir_demux, exist_ok=True)
    os.makedirs(output_dir_demux_fail, exist_ok=True)
    os.makedirs(output_dir_cutadapt, exist_ok=True)

    # demultiplex with cutadapt
    proc_args = [
        "cutadapt",
        "-e",
        "1",
        "--no-indels",
        "-g",
        f"^file:{barcode_fastq}",
        "-o",
        os.path.join(output_dir_demux, "{name}.fq.gz"),
    ]
    if os.path.isdir(fastq_path):
        proc_args.append("-")
        print_command(proc_args)
        cutadapt_demux_proc = subprocess.Popen(
            proc_args,
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        cat_fastq_se(
            fastq_path,
            output_fp=cutadapt_demux_proc.stdin,
        )
        cutadapt_demux_proc.stdin.close()
        cutadapt_demux_proc.wait()
    else:
        proc_args.append(fastq_path)
        print_command(proc_args)
        subprocess.run(proc_args)

    # no need for renaming
    # rename_files_with_mmv(output_dir_demux, rename_pattern)

    subprocess.run(
        [
            "mv",
            os.path.join(output_dir_demux, "unknown.fq.gz"),
            os.path.join(output_dir, "demux_failed"),
        ]
    )

    a, _ = get_primer_set(primer_set)
    cutadapt_trim_proc_args = [
        "cutadapt",
        # "-a" if "..." in a else "-g",
        "-g",
        a,
        "-o",
        output_fastq,
        "--untrimmed-output",
        os.path.join(output_dir_cutadapt, "untrimmed.fq.gz"),
        # "--too-short-output",
        # os.path.join(output_dir_cutadapt, "too_short.fq.gz"),
        "--cores",
        "8",
        "-",
    ]
    if first_k is not None:
        cutadapt_trim_proc_args.extend(["-l", str(first_k)])
    if min_length is not None:
        cutadapt_trim_proc_args.extend(
            [
                "--minimum-length",
                str(min_length),
                "--too-short-output",
                os.path.join(output_dir_cutadapt, "too_short_1.fq.gz"),
            ]
        )
    print_command(cutadapt_trim_proc_args)
    cutadapt_trim_proc = subprocess.Popen(
        cutadapt_trim_proc_args,
        stdin=subprocess.PIPE,
        # silence cutadapt's output
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    cat_fastq_se(
        output_dir_demux,
        output_fp=cutadapt_trim_proc.stdin,
        _have_sample_name=True,
    )
    cutadapt_trim_proc.stdin.close()
    cutadapt_trim_proc.wait()

    shutil.rmtree(output_dir_demux)


def cutadapt_demux_merge_trim_pe(
    fastq_path: str, output_dir: str, primer_set: str, barcode_fastq: str
) -> None:
    """Demultiplex and merge paired-end reads using cutadapt."""
    if not os.path.isfile(barcode_fastq):
        raise ValueError(f"{barcode_fastq} does not exist.")
    output_dir_cutadapt = os.path.join(output_dir, "cutadapt")
    output_dir_demux = os.path.join(output_dir, "demux")
    output_dir_demux_fail = os.path.join(output_dir, "demux_failed")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir_demux, exist_ok=True)
    os.makedirs(output_dir_demux_fail, exist_ok=True)
    os.makedirs(output_dir_cutadapt, exist_ok=True)

    # demultiplex with cutadapt
    proc_args = [
        "cutadapt",
        "-e",
        "1",
        "--no-indels",
        "-g",
        f"^file:{barcode_fastq}",
        "-o",
        os.path.join(output_dir_demux, "{name}_R1.fq.gz"),
        "-p",
        os.path.join(output_dir_demux, "{name}_R2.fq.gz"),
    ]
    if os.path.isdir(fastq_path):
        proc_args.extend(["--interleaved", "-"])
        print_command(proc_args)
        cutadapt_demux_proc = subprocess.Popen(
            proc_args,
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        cat_fastq(
            fastq_path,
            output_fp_r1=cutadapt_demux_proc.stdin,
            output_fp_r2=cutadapt_demux_proc.stdin,
        )
        cutadapt_demux_proc.stdin.close()
        cutadapt_demux_proc.wait()
    else:
        print_command(proc_args)
        proc_args.extend([fastq_path])
        subprocess.run(proc_args)

    # no need for renaming
    # rename_files_with_mmv(output_dir_demux, rename_pattern)
    subprocess.run(
        [
            "mv",
            os.path.join(output_dir_demux, "unknown_R1.fq.gz"),
            os.path.join(output_dir_demux, "unknown_R2.fq.gz"),
            os.path.join(output_dir, "demux_failed"),
        ]
    )
    a, A = get_primer_set(primer_set)
    proc_args = [
        "cutadapt",
        "-e",
        "0.15",
        "-a",
        a,
        "-A",
        A,
        "--pair-filter",
        "first",  # discard pair if first read is not trimmed
        "-o",
        os.path.join(output_dir, "merged_R1.fq.gz"),
        "-p",
        os.path.join(output_dir, "merged_R2.fq.gz"),
        "--untrimmed-output",
        os.path.join(output_dir_cutadapt, "untrimmed_1.fq.gz"),
        "--untrimmed-paired-output",
        os.path.join(output_dir_cutadapt, "untrimmed_2.fq.gz"),
        "--cores",
        "4",
        "--interleaved",
        "-",
    ]
    print_command(proc_args)
    cutadapt_trim_proc = subprocess.Popen(
        proc_args,
        stdin=subprocess.PIPE,
        # silence cutadapt's output
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    cat_fastq(
        output_dir_demux,
        output_fp_r1=cutadapt_trim_proc.stdin,
        output_fp_r2=cutadapt_trim_proc.stdin,
        _have_sample_name=True,
    )
    cutadapt_trim_proc.stdin.close()
    cutadapt_trim_proc.wait()

    shutil.rmtree(output_dir_demux)


def cutadapt_merge_trim_se(
    fastq_path: str,
    output_fastq: str,
    primer_set: str,
    first_k: int | None = None,
    min_length: int | None = None,
    _r2: bool = False,
) -> None:
    output_dir, output_f = os.path.split(output_fastq)
    os.makedirs(output_dir, exist_ok=True)
    a, A = get_primer_set(primer_set)
    proc_args = [
        "cutadapt",
        "-g" if primer_set.endswith("_5") else "-a",
        a,
        "-o",
        output_fastq,
        "--cores",
        "8",
    ]
    if first_k is not None:
        proc_args.extend(["-l", str(first_k)])
    if min_length is not None:
        output_dir_cutadapt = os.path.join(output_dir, "cutadapt")
        os.makedirs(output_dir_cutadapt, exist_ok=True)
        proc_args.extend(
            [
                "--minimum-length",
                str(min_length),
                "--too-short-output",
                os.path.join(output_dir_cutadapt, "too_short_1.fq.gz"),
            ]
        )
    if os.path.isdir(fastq_path):
        proc_args.append("-")
        cutadapt_trim_proc = subprocess.Popen(
            proc_args,
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        cat_fastq_se(fastq_path, output_fp=cutadapt_trim_proc.stdin, _r2=_r2)
        cutadapt_trim_proc.stdin.close()
        cutadapt_trim_proc.wait()
    else:
        proc_args.append(fastq_path)
        subprocess.run(proc_args)


def extract_r2_16s_v3v4(fastq_dir: str, output_fastq: str) -> None:
    length = 150
    output_dir, output_f = os.path.split(output_fastq)
    os.makedirs(output_dir, exist_ok=True)
    log_dir = os.path.join(output_dir, "fastp")
    args = (
        f"fastp -w 8 "
        f"--length_required {length} "
        f"--cut_right "
        f"--stdin "
        f"--stdout "
        f"--json {log_dir}/report.json "
        f"--html {log_dir}/report.html "
        f"| seqkit seq -rp "
        f"| cutadapt -o {output_fastq} -l {length} -"
    )
    print_command(args)
    preprocess_proc = subprocess.Popen(
        args,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    cat_fastq_se(fastq_dir, output_fp=preprocess_proc.stdin, _r2=True)
    preprocess_proc.stdin.close()
    preprocess_proc.wait()


def merge(fastq_dir: str, output_dir: str) -> None:
    """Take output directory of Illumina BCL Convert, merge the fastq files and rename
    each read to `@<sample_name> <read_index> <read 1 or read 2>`. Output are written to
    `output_dir/merged_R1.fastq.gz` and `output_dir/merged_R2.fastq.gz`.

    Note:
        The name of the fastq files are supposed to follow the naming convention of
            `<sample_name>_S<sample_number>_L001_R[1/2]_001.fastq.gz`

    Concatenate files by using subprocess to run something like this:

    for fq1 in "$RAW_DATA_DIR"/*_1.fq; do
        fq2="${fq1%_1.fq}_2.fq"
        sample_name=$(basename "$fq1" "_1.fq")

        awk -v sn="$sample_name" 'BEGIN {
            read_index=1
        }
        {
            if(NR%4==1) {
                print "@sample=" sn " " read_index " 1"
                read_index++
            } else {
                print $0
            }
        }' "$fq1" | gzip >> "$TEMP_DIR/all_samples_1.fq.gz"

        awk -v sn="$sample_name" 'BEGIN {
            read_index=1
        }
        {
            if(NR%4==1) {
                print "@sample=" sn " " read_index " 2"
                read_index++
            } else {
                print $0
            }
        }' "$fq2" | gzip >> "$TEMP_DIR/all_samples_2.fq.gz"

    """
    r1_files = glob.glob(os.path.join(fastq_dir, "*_L001_R1_001.fastq.gz"))
    r2_files = glob.glob(os.path.join(fastq_dir, "*_L001_R2_001.fastq.gz"))
    r1_sample_names = [os.path.basename(f).split("_")[0] for f in r1_files]
    r2_sample_names = [os.path.basename(f).split("_")[0] for f in r2_files]
    # check pairing and uniqueness
    if r1_sample_names != r2_sample_names:
        raise ValueError("Sample names in R1 and R2 do not match.")
    if len(set(r1_sample_names)) != len(r1_sample_names):
        raise ValueError("Sample names are not unique.")

    out_r1 = os.path.join(output_dir, "merged_R1.fastq.gz")
    out_r2 = os.path.join(output_dir, "merged_R2.fastq.gz")
    # remove old files if exist
    if os.path.isfile(out_r1):
        os.remove(out_r1)
    if os.path.isfile(out_r2):
        os.remove(out_r2)
    os.makedirs(output_dir, exist_ok=True)

    # merge
    for r1_file, r2_file, sample_name in zip(r1_files, r2_files, r1_sample_names):
        # merge R1
        for read_num, fq_in, fq_out in zip(
            [1, 2], [r1_file, r2_file], [out_r1, out_r2]
        ):
            process = [
                "awk",
                "-v",
                f"sn={sample_name}",
                f"""'BEGIN {{ read_index=1 }} {{ if(NR%4==1) {{ print "@sample=" sn " " read_index " {read_num}"; read_index++ }} else {{ print $0 }} }}'""",
                fq_in,
                "|",
                "gzip",
                ">>",
                fq_out,
            ]
            subprocess.run(" ".join(process), shell=True)


def _merge_fastq_gz(r1: str, r2: str, output: str) -> None:
    """
    Merges two paired-end gzipped FASTQ files by concatenating the sequence and quality lines.

    Args:
        r1 (str): Path to the first gzipped FASTQ file (R1).
        r2 (str): Path to the second gzipped FASTQ file (R2).
        output (str): Path to the output merged gzipped FASTQ file.
    """
    if not r1.endswith(".gz") or not r2.endswith(".gz") or not output.endswith(".gz"):
        raise ValueError("Input and output files must be gzipped FASTQ files.")
    cmd = f"""
paste -d'\\n' \
<(zcat {r1} | awk 'NR%4==1') \
<(paste -d '' <(zcat {r1} | awk 'NR%4==2') <(zcat {r2} | awk 'NR%4==2')) \
<(zcat {r1} | awk 'NR%4==3') \
<(paste -d '' <(zcat {r1} | awk 'NR%4==0') <(zcat {r2} | awk 'NR%4==0')) \
| gzip -c > {output}"""
    print_command(cmd)

    subprocess.run(cmd, shell=True, executable="/bin/bash", check=True)


#     cmd = f"""
# paste <(zcat {r1} | paste - - - -) <(zcat {r2} | paste - - - -) |
# awk -v OFS="\\n" '{{print $1, $2$6, $3, $4$8}}' | gzip > {output}"""

#     try:
#         subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
#     except subprocess.CalledProcessError as e:
#         raise RuntimeError(f"Error merging FASTQ files: {e}")


def rename_output_files(output_dir, pair):
    """Rename the output files generated by trim_galore."""
    r1_old = os.path.join(output_dir, pair[0].replace("_R1.fastq", "_R1_val_1.fq"))
    r1_new = os.path.join(output_dir, pair[0])
    if not os.path.isfile(r1_old):
        print(
            f"ERROR: {r1_old} does not exist. This is not expected and should be investigated."
        )
    else:
        os.rename(r1_old, r1_new)

    r2_old = os.path.join(output_dir, pair[1].replace("_R2.fastq", "_R2_val_2.fq"))
    r2_new = os.path.join(output_dir, pair[1])
    if not os.path.isfile(r2_old):
        print(
            f"ERROR: {r2_old} does not exist. This is not expected and should be investigated."
        )
    else:
        os.rename(r2_old, r2_new)


def resolve_input_path(input_path: str) -> str:
    """Resolve the input path relative to the 'data' directory if not an existing file.

    Args:
        input_path: A path string or None. If not None and not a file, assume it is relative to the 'data' directory.

    Returns:
        The resolved absolute path or None.
    """
    if not os.path.isfile(input_path):
        return os.path.join(os.path.dirname(__file__), "data", input_path)
    return input_path


# if __name__ == "__main__":
def main():
    parser = argparse.ArgumentParser(
        description="Process and rename FASTQ files using trim_galore."
    )
    parser.add_argument(
        "-i",
        "--input_dir",
        type=str,
        help="Input directory containing FASTQ files",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output FASTQ file (could be gzipped) path or directory",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--primer_set",
        type=str,
        choices=["its", "16s", "maps_round1", "maps_round2"],
        help="Primer set used",
    )
    parser.add_argument(
        "-fb",
        "--barcode_fwd",
        type=str,
        default=None,
        help="Barcode fasta file for forward reads",
    )
    parser.add_argument(
        "-rb",
        "--barcode_rev",
        type=str,
        default=None,
        help="Barcode fasta file for reverse reads",
    )
    parser.add_argument(
        "-pt",
        "--pattern",
        type=str,
        help="Pattern file for mmv, used to rename the files",
    )
    parser.add_argument(
        "-m",
        "--mode",
        type=str,
        default="simple",
        choices=[
            "simple",
            "pseudo_merge",
            "r1",
            "r1_demux_trim",
            "isolate_150",
            "isolate_150_early_stop",
            "maps_round0",
            "maps_round1",
            "maps_round2",
            "maps_round3",
            "maps_rand_hex_test",
        ],
        help="Processing mode",
    )
    parser.add_argument(
        "-k", "--first_k", type=int, default=None, help="The first k bases to keep"
    )
    parser.add_argument(
        "-l", "--min_length", type=int, default=None, help="Minimum length to keep"
    )
    parser.add_argument(
        "-k2", "--first_k2", type=int, default=None, help="The first k bases to keep"
    )
    parser.add_argument(
        "-l2", "--min_length2", type=int, default=None, help="Minimum length to keep"
    )

    args = parser.parse_args()

    if args.mode == "simple":
        simple_preprocess(args.input_dir, args.output)
    elif args.mode == "pseudo_merge":
        # take first_k1 for read1, and first_k2 for read2, then merge them.
        pseudo_merge(
            args.input_dir,
            args.output,
            threads=8,
            min_l1=args.min_length,
            first_k1=args.first_k,
            min_l2=args.min_length2,
            first_k2=args.first_k2,
        )
    elif args.mode == "isolate_150":
        # for barcode_fwd, barcode_rev and pattern, first look for them as if they are
        # normal file path, if not exist, then look for them in script_dir/../../data.
        barcode_fwd = resolve_input_path(args.barcode_fwd)
        barcode_rev = resolve_input_path(args.barcode_rev)
        pattern = resolve_input_path(args.pattern)
        isolate_150_preprocess(
            args.input_dir,
            barcode_fwd,
            barcode_rev,
            pattern,
            args.output,
            primer_set=args.primer_set,
            first_k=args.first_k,
            min_length=args.min_length,
        )
    elif args.mode == "isolate_150_early_stop":
        isolate_150_preprocess(
            args.input_dir,
            args.barcode_fwd,
            args.barcode_rev,
            args.pattern,
            args.output,
            primer_set=args.primer_set,
            first_k=args.first_k,
            min_length=args.min_length,
            early_stop=True,
        )
    elif args.mode == "r1":
        cutadapt_merge_trim_se(
            args.input_dir,
            args.output,
            primer_set=args.primer_set,
            first_k=args.first_k,
            min_length=args.min_length,
        )
    elif args.mode == "r1_demux_trim":
        barcode_fwd = resolve_input_path(args.barcode_fwd)
        cutadapt_demux_merge_trim_se(
            args.input_dir,
            args.output,
            primer_set=args.primer_set,
            barcode_fastq=barcode_fwd,
            first_k=args.first_k,
            min_length=args.min_length,
        )
    elif args.mode == "maps_round1":
        cutadapt_demux_merge_trim_se(
            args.input_dir,
            args.output,
            barcode_fastq=args.barcode_fwd,
            primer_set="maps_1",
        )
    elif args.mode == "maps_round2":
        cutadapt_demux_merge_trim_se(
            args.input_dir,
            args.output,
            barcode_fastq=args.barcode_fwd,
            primer_set="maps_2",
        )
    elif args.mode == "maps_round3":
        cutadapt_demux_merge_trim_se(
            args.input_dir,
            args.output,
            barcode_fastq=args.barcode_fwd,
            primer_set="maps_3",
        )
    elif args.mode == "maps_rand_hex_test":
        cutadapt_demux_merge_trim_pe(
            args.input_dir,
            args.output,
            barcode_fastq=args.barcode_fwd,
            primer_set="maps_rand_hex_test",
        )
