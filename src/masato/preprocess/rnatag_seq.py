import glob
import os
import shutil
import subprocess

from ..trim import rename_files_with_mmv, set_nofile_limit
from ..utils import print_command, cat_fastq, smart_open
from .basic import combine_trim_clip_pe


def rngtagseq_150_preprocess(
    fastq_dir: str,
    barcode_fasta: str,
    rename_pattern: str,
    output_fastq_r1: str,
    output_fastq_r2: str,
    cores: int = 16,
) -> None:
    output_dir, output_f = os.path.split(output_fastq_r1)
    output_dir_cutadapt = os.path.join(output_dir, "cutadapt")
    output_dir_demux = os.path.join(output_dir, "demux")
    output_dir_demux_fail = os.path.join(output_dir, "demux_failed")
    os.makedirs(output_dir_cutadapt, exist_ok=True)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir_demux, exist_ok=True)
    os.makedirs(output_dir_demux_fail, exist_ok=True)

    fastp_json = os.path.join(output_dir_cutadapt, "fastp.json")
    fastp_html = os.path.join(output_dir_cutadapt, "fastp.html")
    trim_args = [
        "fastp",
        "--stdin",
        "--interleaved_in",
        "--stdout",
        "--thread",
        str(cores),
        "--cut_tail",
        "--correction",
        "--html",
        fastp_html,
        "--json",
        fastp_json,
    ]

    cutadapt_args = [
        "cutadapt",
        "-e",
        "0.15",
        "-a",
        "AGATCGGAAGAGCACACGTC;min_overlap=15",
        "-A",
        "AGATCGGAAGAGCGTCGTGT;min_overlap=15",
        "-j",
        str(cores),
        "--interleaved",
        "-",
    ]
    demux_args = [
        "cutadapt",
        "-e",
        "1",
        "-G",
        f"^file:{barcode_fasta}",
        "-j",
        str(cores),
        "--interleaved",
        "-o",
        os.path.join(output_dir_demux, "{name1}-{name2}_R1.fq.gz"),
        "-p",
        os.path.join(output_dir_demux, "{name1}-{name2}_R2.fq.gz"),
        "-",
    ]

    print_command(trim_args)
    trim_proc = subprocess.Popen(
        trim_args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    print_command(cutadapt_args)
    cutadapt_proc = subprocess.Popen(
        cutadapt_args,
        stdin=trim_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    print_command(demux_args)
    demux_proc = subprocess.Popen(
        demux_args,
        stdin=cutadapt_proc.stdout,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    cat_fastq(fastq_dir, output_fp_r1=trim_proc.stdin, output_fp_r2=trim_proc.stdin)
    # ensure all processes are finished
    trim_proc.stdin.close()
    trim_proc.wait()
    # cutadapt_proc.stdin.close()
    cutadapt_proc.wait()
    # demux_proc.stdin.close()
    demux_proc.wait()

    rename_files_with_mmv(output_dir_demux, rename_pattern)

    # move all remaining files (whose name1 or name2 is "unknown") to output_dir/demux_failed
    subprocess.run(
        [
            "mv",
            os.path.join(output_dir_demux, "unknown-unknown_R1.fq.gz"),
            os.path.join(output_dir_demux, "unknown-unknown_R2.fq.gz"),
            os.path.join(output_dir, "demux_failed"),
        ]
    )

    awk_code = r"""{
    if (NR % 4 == 1) {
        split($0, parts, " ")
        sub(/^@/, "", parts[1])
        suffix = substr(parts[1], 8)  # skip first 7 characters
        new_id = parts[4] "_" suffix
        print "@" new_id
    } else {
        print
    }
}"""
    awk_args = ["awk", awk_code]
    print_command(awk_args)
    awk1_proc = subprocess.Popen(
        awk_args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        # stderr=subprocess.DEVNULL,
    )
    seqkit1_args = ["seqkit", "seq", "-o", output_fastq_r1]
    print_command(seqkit1_args)
    seqkit1_proc = subprocess.Popen(
        seqkit1_args,
        stdin=awk1_proc.stdout,
        stdout=subprocess.DEVNULL,
        # stderr=subprocess.DEVNULL,
    )
    print_command(awk_args)
    awk2_proc = subprocess.Popen(
        awk_args,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL,
    )
    seqkit2_args = ["seqkit", "seq", "-o", output_fastq_r2]
    print_command(seqkit2_args)
    seqkit2_proc = subprocess.Popen(
        seqkit2_args,
        stdin=awk2_proc.stdout,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    cat_fastq(
        output_dir_demux,
        output_fp_r1=awk1_proc.stdin,
        output_fp_r2=awk2_proc.stdin,
        _have_sample_name=True,
    )
