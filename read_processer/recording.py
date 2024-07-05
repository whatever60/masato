import os
import subprocess
import gzip

from Bio import SeqIO

import sys

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from trim import get_primer_set, get_rc
from utils import print_command, cat_fastq, find_paired_end_files


def run_fastp(sample: str, read1: str, read2: str, output_dir: str, cpus: int):
    [
        os.makedirs(os.path.join(output_dir, d), exist_ok=True)
        for d in ["merged", "unpaired", "unmerged", "failed", "log"]
    ]

    merge_cmd = [
        "fastp",
        "-i",
        read1,
        "-I",
        read2,
        "-w",
        str(cpus),
        "--merge",
        "--merged_out",
        f"{output_dir}/merged/{sample}.fq.gz",
        "--unpaired1",
        f"{output_dir}/unpaired/{sample}.1.fq.gz",
        "--unpaired2",
        f"{output_dir}/unpaired/{sample}.2.fq.gz",
        "--out1",
        f"{output_dir}/unmerged/{sample}.1.fq.gz",
        "--out2",
        f"{output_dir}/unmerged/{sample}.2.fq.gz",
        "--failed_out",
        f"{output_dir}/failed/{sample}.fq.gz",
        "--html",
        f"{output_dir}/log/{sample}.html",
        "--json",
        f"{output_dir}/log/{sample}.json",
    ]

    with open(f"{output_dir}/log/{sample}.merge.err", "w") as err_file:
        print_command(merge_cmd)
        subprocess.run(merge_cmd, stderr=err_file)


def run_cutadapt(
    input_file: str,
    output_file: str,
    *,
    untrimmed_output: str = None,
    log_file: str = None,
    adapter_5: str = None,
    adapter_3: str = None,
    cores: int,
    e: int | float = 0.05,
    rc: bool = False,
):
    # makedirs for all output locations
    [
        os.makedirs(os.path.dirname(f), exist_ok=True)
        for f in [output_file, untrimmed_output, log_file]
        if f is not None
    ]

    def add_adapter(cmd: list):
        if adapter_5 is not None:
            cmd.extend(["-g", adapter_5])
        if adapter_3 is not None:
            cmd.extend(["-a", adapter_3])

    if rc:
        seqkit_cmd_1 = [
            "seqkit",
            "seq",
            "-r",
            "-p",
            "-t",
            "DNA",
            "--validate-seq",
            input_file,
        ]
        cutadapt_cmd = [
            "cutadapt",
            "-e",
            str(e),
            "--cores",
            str(cores),
        ]
        add_adapter(cutadapt_cmd)
        if untrimmed_output is not None:
            cutadapt_cmd.extend(["--untrimmed-output", untrimmed_output])
        cutadapt_cmd.append("-")
        seqkit_cmd_2 = [
            "seqkit",
            "seq",
            "-r",
            "-p",
            "-t",
            "DNA",
            "--validate-seq",
            "-",
            "-o",
            output_file,
        ]

        with open(log_file, "w") as log:
            print_command(cutadapt_cmd)
            p1 = subprocess.Popen(seqkit_cmd_1, stdout=subprocess.PIPE)
            print_command(cutadapt_cmd)
            p2 = subprocess.Popen(
                cutadapt_cmd,
                stdin=p1.stdout,
                stdout=subprocess.PIPE,
                stderr=log,
            )
            print_command(seqkit_cmd_2)
            p3 = subprocess.Popen(seqkit_cmd_2, stdin=p2.stdout, stdout=subprocess.PIPE)
            p3.communicate()

    else:
        cutadapt_cmd = [
            "cutadapt",
            "-e",
            str(e),
            "-o",
            output_file,
            "--cores",
            str(cores),
        ]
        add_adapter(cutadapt_cmd)

        if untrimmed_output:
            cutadapt_cmd.extend(["--untrimmed-output", untrimmed_output])

        cutadapt_cmd.append(input_file)

        with open(log_file, "w") as log:
            print_command(cutadapt_cmd)
            subprocess.run(cutadapt_cmd, stdout=log)


def makedirs(*directories):
    [os.makedirs(d, exist_ok=True) for d in directories]


def process_one_sample(sample: str, read1: str, read2: str, output_dir: str, cpus: int):
    merge_dir = os.path.join(output_dir, "merge")
    os.makedirs(merge_dir, exist_ok=True)


    # Merge pairs and remove adapters
    output_mp = f"{merge_dir}/merged/{sample}.fq.gz"
    output_rm_adapters = f"{merge_dir}/cleaned/{sample}.fq.gz"
    run_fastp(sample, read1, read2, merge_dir, cpus)
    run_cutadapt(
        input_file=output_mp,
        output_file=output_rm_adapters,
        log_file=f"{merge_dir}/log/{sample}.clean.out",
        adapter_3=get_primer_set("recording_adapter"),
        cores=cpus,
    )

    # 5' UMI
    umi5_dir = os.path.join(output_dir, "umi5")
    makedirs(umi5_dir)
    output_umi5_rest = f"{umi5_dir}/rest/{sample}.fq.gz"
    output_umi5_umi = f"{umi5_dir}/umi/{sample}.fq.gz"
    run_cutadapt(
        input_file=output_rm_adapters,
        output_file=output_umi5_umi,
        untrimmed_output=f"{umi5_dir}/umi_untrimmed/{sample}.fq.gz",
        log_file=f"{umi5_dir}/log/{sample}.umi.out",
        adapter_3=get_primer_set("recording_leader_dr"),
        cores=cpus,
    )
    run_cutadapt(
        input_file=output_rm_adapters,
        output_file=output_umi5_rest,
        untrimmed_output=f"{umi5_dir}/rest_untrimmed/{sample}.fq.gz",
        log_file=f"{umi5_dir}/log/{sample}.rest.out",
        adapter_5=get_primer_set("recording_leader_dr"),
        cores=cpus,
    )

    # 3' UMI
    umi3_dir = os.path.join(output_dir, "umi3")
    makedirs(umi3_dir)
    output_umi3_rest = f"{umi3_dir}/rest/{sample}.fq.gz"
    output_umi3_umi = f"{umi3_dir}/umi/{sample}.fq.gz"
    run_cutadapt(
        input_file=output_umi5_rest,
        output_file=output_umi3_umi,
        untrimmed_output=f"{umi3_dir}/umi_untrimmed/{sample}.fq.gz",
        log_file=f"{umi3_dir}/log/{sample}.umi.out",
        adapter_3=get_primer_set("recording_spacer0") + ";rightmost",
        cores=cpus,
    )
    run_cutadapt(
        input_file=output_umi5_rest,
        output_file=output_umi3_rest,
        untrimmed_output=f"{umi3_dir}/rest_untrimmed/{sample}.fq.gz",
        log_file=f"{umi3_dir}/log/{sample}.rest.out",
        adapter_5=get_rc(get_primer_set("recording_spacer0")),
        cores=cpus,
        rc=True,
    )

    # Spacer rounds
    round_num = 1
    current_input = output_umi3_rest
    while True:
        spacer_dir = os.path.join(output_dir, f"spacer{round_num}")
        makedirs(spacer_dir)
        output_spacer_rest = f"{spacer_dir}/rest/{sample}.fq.gz"
        output_spacer_spacer = f"{spacer_dir}/spacer/{sample}.fq.gz"
        run_cutadapt(
            input_file=current_input,
            output_file=output_spacer_spacer,
            untrimmed_output=f"{spacer_dir}/spacer_untrimmed/{sample}.fq.gz",
            log_file=f"{spacer_dir}/log/{sample}.spacer.out",
            adapter_3=get_primer_set("recording_dr"),
            cores=cpus,
        )
        run_cutadapt(
            input_file=current_input,
            output_file=output_spacer_rest,
            untrimmed_output=f"{spacer_dir}/rest_untrimmed/{sample}.fq.gz",
            log_file=f"{spacer_dir}/log/{sample}.rest.out",
            adapter_5=get_primer_set("recording_dr"),
            cores=cpus,
        )
        if all_sequences_empty(output_spacer_rest):
            break

        current_input = output_spacer_rest
        round_num += 1


def all_sequences_empty(fastq_file: str) -> bool:
    """Check if all sequences in a FASTQ file are empty, or if the file itself is empty."""
    with gzip.open(fastq_file, "rt") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            if len(record.seq) > 0:
                return False
    return True


if __name__ == "__main__":
    process_one_sample(
        "20240628RM29",
        read1="/mnt/c/aws_data/20240701_yuanyuan_recording/fastq_raw/20240628RM29_S29_L001_R1_001.fastq.gz",
        read2="/mnt/c/aws_data/20240701_yuanyuan_recording/fastq_raw/20240628RM29_S29_L001_R2_001.fastq.gz",
        output_dir="/mnt/c/aws_data/20240701_yuanyuan_recording/read_process_test",
        cpus=4,
    )
