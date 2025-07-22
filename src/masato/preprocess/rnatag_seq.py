import os
import subprocess
import tempfile

from Bio import SeqIO
from tqdm.auto import tqdm

from ..trim import TRUSEQ_READ1, TRUSEQ_READ2, get_rc
from ..utils import print_command, find_paired_end_files


def rnatagseq_150_preprocess(
    input_: str | list[str],
    barcode_fasta: str,
    output_dir: str,
    cores: int = 16,
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    truseq_adapter_config = "min_overlap=5;e=0.15"

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta") as barcode_file:
        for record in SeqIO.parse(barcode_fasta, "fasta"):
            # write name
            seq = f"^{record.seq};e=1...r1={get_rc(TRUSEQ_READ1)};{truseq_adapter_config};optional"
            barcode_file.write(f">{record.id}\n{seq}\n")
        barcode_file.flush()

        prog_bar = tqdm(find_paired_end_files(input_))
        for r1_path, r2_path, sample_name in prog_bar:
            prog_bar.set_description(f"Processing {sample_name}")
            output_dir_sample = os.path.join(output_dir, sample_name)
            os.makedirs(output_dir_sample, exist_ok=True)
            output_r1 = os.path.join(output_dir_sample, "merged_1.fq.gz")
            output_r2 = os.path.join(output_dir_sample, "merged_2.fq.gz")
            fastp_json = os.path.join(output_dir_sample, "fastp.json")
            fastp_html = os.path.join(output_dir_sample, "fastp.html")

            trim_args = [
                "fastp",
                "--in1",
                r1_path,
                "--in2",
                r2_path,
                "--stdout",
                "--thread",
                str(cores),
                "--cut_tail",
                "--correction",
                "--disable_length_filtering",  # default is 15 and will filter out 8bp read2.
                "--html",
                fastp_html,
                "--json",
                fastp_json,
            ]

            # NOTE about cutadapt 3' adapter removal and 5' barcode demultiplexing at the same time:
            # cutadapt only trims each read n times (--times) and the default is 1. So by default,
            # only 5' OR 3' adapter will be removed, not both. However, cutadapt does go through
            # all adapters and choose the best match for trimming. So to trim both 5' and 3' adapters,
            # we can consider the following options:
            # 1. Use linked adatper. It would require an awkward fasta file. But I am adopting
            # it here while requesting a new feature at https://github.com/marcelm/cutadapt/issues/851.
            # 2. Run cutadapt twice. It's easier to remove reads with bad barcodes by discarding untrimmed.
            # 3. Specify both -A and -G and --times 2. This is only approximately correct.

            # NOTE about cutadapt anchored adapters:
            # Anchored adapters are required to fully match. So min_overlap is not allowed or
            # ignored. Therefore, to accommodate sequencing configs where the second read is only
            # as long as the barcode (i.e., linker is not sequences), we need to allow indels,
            # i.e., we cannot specify `noindels`, even that will be more rigid and faster.

            # NOTE: Create a temp file from a legit barcode fasta file by doing
            # f"^<seq>;e=1...r1=<get_rc(truseq_read1)>;min_overlap=5;e=0.15;optional"
            cutadapt_args = [
                "cutadapt",
                "-a",
                f"r2={TRUSEQ_READ2};{truseq_adapter_config}",
                "-A",
                f"file:{barcode_file.name}",
                "--rename",
                "{id}_{r2.adapter_name} {comment}",
                "-j",
                str(cores),
                "--interleaved",
                "-o",
                output_r1,
                "-p",
                output_r2,
                "-",
            ]

            print_command(trim_args)
            trim_proc = subprocess.Popen(
                trim_args,
                # stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.DEVNULL,
            )
            print_command(cutadapt_args)
            cutadapt_proc = subprocess.Popen(
                cutadapt_args,
                stdin=trim_proc.stdout,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            trim_proc.wait()
            cutadapt_proc.wait()
