import subprocess


def map_se(
    reads_fastq: str,
    reference_genome: str,
    *,
    bam_file: str,
    bwa_log: str = "/dev/null",
    mapped_reads: str,
    unmapped_reads: str,
    args: str = "default",
    num_threads: int = 12,
):
    """
    Aligns reads to a reference genome, sorts, indexes the BAM file, and extracts mapped and unmapped reads.

    Args:
        reference_genome (str): Path to the reference genome file.
        reads_fastq (str): Path to the input FASTQ file.
        mapped_reads (str): Path to the output file for mapped reads.
        unmapped_reads (str): Path to the output file for unmapped reads.
        bam_file (str): Path to the sorted BAM file.
        bwa_log (str, optional): Path to the file where stderr of bwa will be logged. Defaults to None.
    """
    if args == "default":
        args = ""
    elif args == "relaxed":
        args = "-k 9 -B 2 -O 4 -E 2 -T 27 "
        # B: mismatch penalty, O: gap open penalty, E: gap extension penalty, T: minimum threshold
    try:
        # Step 2: Align the reads and pipe directly to BAM, then sort
        if not bwa_log:
            bwa_log = "/dev/null"
        align_cmd = (
            f"bwa-mem2 mem -t {num_threads} {args}{reference_genome} {reads_fastq} "
            f"2> {bwa_log} | samtools view -S -b - | samtools sort -o {bam_file}"
        )
        _ = subprocess.run(align_cmd, shell=True, check=True)

        # Step 3: Index the sorted BAM file
        index_cmd = f"samtools index {bam_file}"
        _ = subprocess.run(index_cmd, shell=True, check=True)

        # Step 4: Extract unmapped reads
        extract_mapped_cmd = (
            f"samtools view -F 4 -b {bam_file} | samtools fastq - | pigz > {mapped_reads}"
        )
        _ = subprocess.run(
            extract_mapped_cmd, shell=True, check=True, stderr=subprocess.DEVNULL
        )

        extract_unmapped_cmd = (
            f"samtools view -f 4 -b {bam_file} | samtools fastq - | pigz > {unmapped_reads}"
        )
        _ = subprocess.run(
            extract_unmapped_cmd, shell=True, check=True, stderr=subprocess.DEVNULL
        )

    except subprocess.CalledProcessError as e:
        print(f"An error occurred while processing reads: {e}")


# Example usage:
# process_reads("path/to/ecrec.fna", "path/to/spacers.fq.gz", "path/to/spacers_self.fa", "path/to/spacers_others.fa", "path/to/sorted_aligned_reads.bam", "path/to/bwa_log.txt")
