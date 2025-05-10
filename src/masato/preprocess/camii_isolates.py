import glob
import os
import shutil
import subprocess

from ..trim import get_primer_set, rename_files_with_mmv, set_nofile_limit
from ..utils import print_command, cat_fastq


def isolate_150_preprocess(
    fastq_dir: str,
    barcode_fwd_fasta: str,
    barcode_rev_fasta: str,
    rename_pattern: str,
    output_fastq: str,
    primer_set: str,
    first_k: int | None = None,
    first_k_r2: int | None = None,
    min_length: int = 0,
    min_length_r2: int = 0,
    early_stop: bool = False,
    quality_trimming: bool = False,
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
        preexec_fn=set_nofile_limit,
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
        ["cutadapt", "-a", a]
        + (["-A", A] if not early_stop else [])  # no trimming for read 2 if early stop
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
        proc_args.extend(["--length", str(first_k)])
    if first_k_r2 is not None:
        proc_args.extend(["-L", str(first_k_r2)])
    print_command(proc_args)
    cutadapt_trim_proc = subprocess.Popen(
        proc_args,
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
        _fastp_trim_proc = subprocess.Popen(
            proc_args,
            stdin=subprocess.PIPE,
            stdout=_cutadapt_trim_proc.stdin,
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

    # copy merged_1 to output_fastq
    shutil.move(os.path.join(output_dir_cutadapt, output_fastq_r1), output_fastq)
    shutil.move(
        os.path.join(output_dir_cutadapt, output_fastq_r2),
        os.path.join(os.path.dirname(output_fastq), output_fastq_r2),
    )
    shutil.rmtree(output_dir_demux)
