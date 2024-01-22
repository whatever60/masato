#!/usr/bin/env python

import glob
import os
import subprocess
import argparse

from Bio.Seq import Seq

from utils import cat_fastq


PRIMER_ITS_5 = "ACCTGCGGARGGATCA"
PRIMER_ITS_7 = "GAGATCCRTTGYTRAAAGTT"
PRIMER_16S_5 = "GTGYCAGCMGCCGCGGTAA"
PRIMER_16S_7 = "GGACTACNVGGGTWTCTAAT"


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
) -> None:
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
    output_dir: str = None,
    amplicon_type: str = "ITS",
):
    """Quality trimming and filtering by fastp, then adapter trimming by cutadapt.
    For read 1, trim RC of i7 adapter from the 3' end; For read 2, trim RC of i5 adapter
        also from the 3' end. Just trim and save report, no need to separate trimmed
        and untrimmed reads or discard reads.
    """
    output_fastq_1 = os.path.join(output_dir, os.path.basename(input_fastq_1))
    output_fastq_2 = os.path.join(output_dir, os.path.basename(input_fastq_2))
    report_report = os.path.join(output_dir, "cutadapt.report.json")
    if not amplicon_type in ["ITS", "16S"]:
        raise ValueError("amplicon_type must be either ITS or 16S.")

    fastp_proc = run_fastp(input_fastq_1, input_fastq_2, output_dir, stdout=True)
    subprocess.run(
        [
            "cutadapt",
            "-a",
            str(
                Seq(
                    PRIMER_ITS_7 if amplicon_type == "ITS" else PRIMER_16S_7
                ).reverse_complement()
            ),
            "-A",
            str(
                Seq(
                    PRIMER_ITS_5 if amplicon_type == "ITS" else PRIMER_16S_5
                ).reverse_complement()
            ),
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
            "100",
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process and rename FASTQ files using trim_galore."
    )
    parser.add_argument(
        "-i", "--input_dir", type=str, help="Input directory containing FASTQ files"
    )
    parser.add_argument(
        "-o",
        "--output_fastq",
        type=str,
        help="Output FASTQ file path, could be gzipped",
    )
    args = parser.parse_args()

    simple_preprocess(args.input_dir, args.output_fastq)
