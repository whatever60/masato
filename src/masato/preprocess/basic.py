import os
import subprocess

from ..trim import get_primer_set
from ..utils import print_command, cat_fastq


def combine_trim_clip_pe(
    input_fastq_dir: str,
    output_fastq: str,
    output_fastq_r2: str,
    primer_set: str | None = None,
    first_k: int | None = None,
    first_k_r2: int | None = None,
    min_length: int = 0,
    min_length_r2: int = 0,
    quality_trimming: bool = False,
    cores: int = 16,
    _have_sample_name: bool = False,
) -> None:
    """
    Combine, quality trim, adapter trim, clip, and length-filter paired-end reads.

    1. Combine using custom function.
    2. (Optional) Quality trim using fastq.
    3. Trim adapter, clip to maximum length and filter by length using cutadapt.
    """
    output_dir, _ = os.path.split(output_fastq)
    os.makedirs(output_dir, exist_ok=True)
    output_dir_cutadapt = os.path.join(output_dir, "cutadapt")
    os.makedirs(output_dir_cutadapt, exist_ok=True)
    if primer_set is not None:
        a, A = get_primer_set(primer_set)
        proc_args = [
            "cutadapt",
            "-e",
            "0.15",
            "-O",
            "16",
            "-g",
            a,
            "-G",
            A,
            "--untrimmed-output",
            os.path.join(output_dir_cutadapt, "untrimmed_1.fq.gz"),
            "--untrimmed-paired-output",
            os.path.join(output_dir_cutadapt, "untrimmed_2.fq.gz"),
        ]
    else:
        proc_args = ["cutadapt"]
    if first_k is not None:
        proc_args += ["--length", str(first_k)]
    if first_k_r2 is not None:
        proc_args += ["-L", str(first_k_r2)]
    proc_args += [
        "--pair-filter",
        "any",
        "--minimum-length",
        f"{str(min_length)}:{str(min_length_r2)}",
        "--too-short-output",
        os.path.join(output_dir_cutadapt, "too_short_1.fq.gz"),
        "--too-short-paired-output",
        os.path.join(output_dir_cutadapt, "too_short_2.fq.gz"),
        "--cores",
        str(cores),
        "--interleaved",
        "-o",
        output_fastq,
        "-p",
        output_fastq_r2,
        "-",
    ]
    print_command(proc_args)
    cutadapt_trim_proc = subprocess.Popen(
        proc_args,
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    if quality_trimming:
        _cutadapt_trim_proc = cutadapt_trim_proc
        fastp_json = os.path.join(output_dir_cutadapt, "fastp.json")
        fastp_html = os.path.join(output_dir_cutadapt, "fastp.html")
        proc_args = [
            "fastp",
            "--stdin",
            "--interleaved_in",
            "--stdout",
            "--cut_right",
            "--correct",
            "--thread",
            "1",  # Threads can only be 1 for consistent output order as input
            "--json",
            fastp_json,
            "--html",
            fastp_html,
        ]
        print_command(proc_args)
        _fastp_trim_proc = subprocess.Popen(
            proc_args,
            stdin=subprocess.PIPE,
            stdout=cutadapt_trim_proc.stdin,
            stderr=subprocess.DEVNULL,
        )
        cutadapt_trim_proc = _fastp_trim_proc

    cat_fastq(
        input_fastq_dir,
        output_fp_r1=cutadapt_trim_proc.stdin,
        output_fp_r2=cutadapt_trim_proc.stdin,
        _have_sample_name=_have_sample_name,
    )
    cutadapt_trim_proc.stdin.close()
    cutadapt_trim_proc.wait()


def combine_trim_merge_pe(
    input_dir: str,
    output_fastq: str,
    min_length: int,
    cores: int,
) -> None:
    """
    Given a directory of fastq files, perform renaming, quality trimming, adapter
    trimming, length filtering, and merging using fastp. Fastp takes interleaved
    input from stdin, to which we write the renamed and interleaved paired-end reads.
    """
    output_dir = os.path.dirname(output_fastq)
    os.makedirs(output_dir, exist_ok=True)
    fastp_output_dir = os.path.join(output_dir, "fastp")
    os.makedirs(fastp_output_dir, exist_ok=True)

    proc_args = [
        "fastp",
        "--stdin",
        "--interleaved_in",
        "--thread",  # Use --thread as requested
        str(cores),
        "--length_required",
        str(min_length),
        "--cut_right",
        "--merge",
        "--merged_out",
        output_fastq,
        # --- Output for non-merged/failed reads ---
        "--unpaired1",
        os.path.join(fastp_output_dir, "unpaired_R1.fastq.gz"),
        "--unpaired2",
        os.path.join(fastp_output_dir, "unpaired_R2.fastq.gz"),
        "--failed_out",
        os.path.join(fastp_output_dir, "failed.fastq.gz"),
        "--out1",
        os.path.join(fastp_output_dir, "unmerged_R1.fastq.gz"),
        "--out2",
        os.path.join(fastp_output_dir, "unmerged_R2.fastq.gz"),
        # --- Reporting ---
        "--html",
        os.path.join(fastp_output_dir, "report.html"),
        "--json",
        os.path.join(fastp_output_dir, "report.json"),
    ]

    fastp_proc = subprocess.Popen(
        proc_args,
        stdin=subprocess.PIPE,
    )

    # Assuming `cat_fastq` is defined elsewhere and handles writing to the process stdin
    cat_fastq(input_dir, output_fp_r1=fastp_proc.stdin, output_fp_r2=fastp_proc.stdin)

    fastp_proc.stdin.close()
    fastp_proc.wait()
