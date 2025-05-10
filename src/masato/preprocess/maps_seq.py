import os
import subprocess
import shutil

from ..trim import get_primer_set, set_nofile_limit
from ..utils import print_command, cat_fastq, cat_fastq_se


def cutadapt_demux_merge_trim_se(
    fastq_path: str,
    output_fastq: str,
    output_fastq_r2: str,
    primer_set: str,
    barcode_fastq: str,
    first_k: int | None = None,
    first_k_r2: int | None = None,
    min_length: int = 0,
    min_length_r2: int = 0,
    quality_trimming: bool = False,
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
        f"^file:{barcode_fastq}",  # There is no -G in this function
        "-o",
        os.path.join(output_dir_demux, "{name}_R1.fq.gz"),
        "-p",
        os.path.join(output_dir_demux, "{name}_R2.fq.gz"),
        "--cores",
        "1",
        "--interleaved",
        "-",
    ]
    if os.path.isdir(fastq_path):
        proc_args.append("-")
        print_command(proc_args)
        cutadapt_demux_proc = subprocess.Popen(
            proc_args,
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            preexec_fn=set_nofile_limit,
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

    a, A = get_primer_set(primer_set)
    cutadapt_trim_proc_args = (
        ["cutadapt", "-g", a]
        + (["-G", A] if A is not None else [])
        + [
            "-e",
            "0.15",
            "-O",
            "16",
            "--pair-filter",
            "any",
            "--minimum-length",
            f"{str(min_length)}:{str(min_length_r2)}",
            "-o",
            output_fastq,
            "-p",
            output_fastq_r2,
            "--untrimmed-output",
            os.path.join(output_dir_cutadapt, "untrimmed_1.fq.gz"),
            "--untrimmed-paired-output",
            os.path.join(output_dir_cutadapt, "untrimmed_2.fq.gz"),
            "--too-short-output",
            os.path.join(output_dir_cutadapt, "too_short_1.fq.gz"),
            "--too-short-paired-output",
            os.path.join(output_dir_cutadapt, "too_short_2.fq.gz"),
            "--cores",
            "8",
            "--interleaved",
            "-",
        ]
    )
    if first_k is not None:
        cutadapt_trim_proc_args.extend(["--length", str(first_k)])
    if first_k_r2 is not None:
        cutadapt_trim_proc_args.extend(["-L", str(first_k_r2)])
    print_command(cutadapt_trim_proc_args)
    cutadapt_trim_proc = subprocess.Popen(
        cutadapt_trim_proc_args,
        stdin=subprocess.PIPE,
        # silence cutadapt's output
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
            "--thread",
            "8",
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
        output_dir_demux,
        output_fp_r1=cutadapt_trim_proc.stdin,
        output_fp_r2=cutadapt_trim_proc.stdin,
        _have_sample_name=True,
    )
    cutadapt_trim_proc.stdin.close()
    cutadapt_trim_proc.wait()

    shutil.rmtree(output_dir_demux)
