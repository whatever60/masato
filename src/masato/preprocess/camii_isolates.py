import glob
import os
import shutil
import subprocess

from ..trim import rename_files_with_mmv, set_nofile_limit
from ..utils import print_command, cat_fastq
from .basic import combine_trim_clip_pe


def isolate_150_preprocess(
    fastq_dir: str,
    barcode_fwd_fasta: str,
    barcode_rev_fasta: str,
    rename_pattern: str,
    output_fastq: str,
    output_fastq_r2: str,
    # args for combine_trim_clip_pe
    primer_set: str,
    first_k: int | None = None,
    first_k_r2: int | None = None,
    min_length: int = 0,
    min_length_r2: int = 0,
    # early_stop: bool = False,
    quality_trimming: bool = False,
    cores: int = 16,
) -> None:
    output_dir, output_f = os.path.split(output_fastq)
    output_dir_demux = os.path.join(output_dir, "demux")
    output_dir_demux_fail = os.path.join(output_dir, "demux_failed")
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(output_dir_demux, exist_ok=True)
    os.makedirs(output_dir_demux_fail, exist_ok=True)

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

    combine_trim_clip_pe(
        input_fastq_dir=output_dir_demux,
        output_fastq=output_fastq,
        output_fastq_r2=output_fastq_r2,
        primer_set=primer_set,
        first_k=first_k,
        first_k_r2=first_k_r2,
        min_length=min_length,
        min_length_r2=min_length_r2,
        quality_trimming=quality_trimming,
        cores=cores,
        _have_sample_name=True,
    )

    shutil.rmtree(output_dir_demux)
