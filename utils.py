import io
import glob
import gzip
import os
import re
import warnings

import pandas as pd
from tqdm.auto import tqdm


def pattern_illumina(file_: str) -> tuple[str, str]:
    """Pattern function for Illumina default output."""
    file_pattern = re.compile(r"^(.+?)_S\d+_L\d{3}_(R[12])_001.f(ast)?q.gz$")
    match = file_pattern.search(os.path.basename(file_))
    if match:
        sample_name, read_type = match.groups()[:2]
        return sample_name, read_type
    return None


def pattern_custom(file_: str) -> tuple[str, str]:
    """Pattern function for custom pattern ending with fq.gz or fastq.gz."""
    file_pattern = re.compile(r"^(.+?)_(R[12]).f(ast)?q.gz$")
    match = file_pattern.search(os.path.basename(file_))
    if match:
        sample_name, read_type = match.groups()[:2]
        return sample_name, read_type
    return None


def find_paired_end_files(directory: str) -> list[tuple[str, str, str]]:
    pattern_functions = [pattern_illumina, pattern_custom]

    # Dictionary to store file pairs
    file_pairs = {}

    # Find all .fastq.gz and .fq.gz files in the directory
    fastq_files = glob.glob(os.path.join(directory, "*.fastq.gz")) + glob.glob(
        os.path.join(directory, "*.fq.gz")
    )

    # Process each file using pattern functions
    for file_ in fastq_files:
        matched = False
        for pattern_function in pattern_functions:
            match = pattern_function(file_)
            if match:
                sample_name, read_type = match
                if read_type not in ["R1", "R2"]:
                    warnings.warn(f"Invalid read type for file: {file_}")
                    continue

                pair_key = (sample_name, read_type)

                # Add the file to the dictionary and check for duplicates
                if pair_key in file_pairs:
                    warnings.warn(f"Duplicate file for {pair_key}: {file_}")
                else:
                    file_pairs[pair_key] = file_
                matched = True
                break

        if not matched:
            warnings.warn(f"File {file_} does not match any known pattern.")

    # Match R1 and R2 pairs and prepare the output
    matched_pairs = []
    processed_samples = set()

    for (sample_name, read_type), file_ in file_pairs.items():
        if sample_name in processed_samples:
            # Skip if we've already processed this sample
            continue

        processed_samples.add(sample_name)
        pair_type = "R2" if read_type == "R1" else "R1"
        pair_key = (sample_name, pair_type)

        if pair_key in file_pairs:
            r1_file = file_ if read_type == "R1" else file_pairs[pair_key]
            r2_file = file_pairs[pair_key] if read_type == "R1" else file_
            matched_pairs.append((r1_file, r2_file, sample_name))
        else:
            warnings.warn(f"Missing pair for file: {file_}")

    # Report the number of matched pairs
    print(f"Number of matched pairs: {len(matched_pairs)}")

    # Order by sample name
    matched_pairs.sort(key=lambda x: x[2])
    return matched_pairs


def cat_fastq(
    directory: str,
    output_fp_r1,
    output_fp_r2=None,
    metadata: str = None,
    _remove_undet: bool = False,
):
    """Process FASTQ files in the given directory, renaming reads,and write the output
    to the specified file pointers. Output fastq will be interleaved if `output_fp_r2`
    is None.

    This function reinvents the wheel implemented in many bioinformatics tools, but I
    am still doing this for customizing the read renaming.

    Args:
        directory: Directory containing FASTQ files.
        output_fp_r1: File pointer to write the R1 output.
        output_fp_r2: File pointer to write the R2 output.
    """
    matched_pairs = find_paired_end_files(directory)
    if output_fp_r2 is None:
        output_fp_r2 = output_fp_r1
    if metadata is not None:
        samples_in_meta = pd.read_table(metadata, index_col=0).index.to_list()
    else:

        class _all_in:
            def __contains__(self, item):
                return True

        samples_in_meta = _all_in()

    if (
        (
            isinstance(output_fp_r1, gzip.GzipFile)
            and isinstance(output_fp_r2, gzip.GzipFile)
        )
        or (
            isinstance(output_fp_r1, io.BufferedWriter)
            and isinstance(output_fp_r2, io.BufferedWriter)
        )
        or (
            isinstance(output_fp_r1, io.BytesIO)
            and isinstance(output_fp_r2, io.BytesIO)
        )
    ):

        def write_line(line1, line2: str):
            output_fp_r1.write(line1.encode())
            output_fp_r2.write(line2.encode())

    elif isinstance(output_fp_r1, io.TextIOBase) and isinstance(
        output_fp_r2, io.TextIOBase
    ):

        def write_line(line1, line2: str):
            output_fp_r1.write(line1)
            output_fp_r2.write(line2)

    else:
        raise ValueError("Output file pointers must be both gzip or both text file.")

    for r1_path, r2_path, sample_name in tqdm(matched_pairs):
        if not sample_name in samples_in_meta:
            continue
        if _remove_undet and sample_name == "Undetermined":
            continue
        with gzip.open(r1_path, "rt") as r1_file, gzip.open(r2_path, "rt") as r2_file:
            paired_read_iter = zip(
                zip(*[r1_file] * 4, strict=True),
                zip(*[r2_file] * 4, strict=True),
                strict=True,
            )

            for read_index, (r1_lines, r2_lines) in enumerate(
                paired_read_iter, start=1
            ):
                # Renaming reads
                write_line(
                    rename_read(r1_lines[0], sample_name, 1, read_index)
                    + "".join(r1_lines[1:]),
                    rename_read(r2_lines[0], sample_name, 2, read_index)
                    + "".join(r2_lines[1:]),
                )


# def reads_in_four(iterable: Iterator[str]) -> Iterator[tuple[str, str, str, str]]:
#     """
#     Yield groups of four lines from the given iterable.

#     Args:
#         iterable: An iterator from which to read lines.

#     Returns:
#         An iterator of tuples, each containing four lines.
#     """
#     iters = [iter(iterable)] * 4
#     return zip(*iters)


def rename_read(
    header_line: str, sample_name: str, read_number: int, read_index: int
) -> str:
    """
    Rename a read header line with the sample name, read number, and read index.

    Args:
        header_line: Original header line from FASTQ file.
        sample_name: Name of the sample.
        read_number: Read number (1 or 2).
        read_index: Index of the read in the sample.

    Returns:
        Renamed header line.
    """
    # Remove '@' and split by space, take first part
    original_header = header_line.strip().split()[0][1:]
    new_header = f"@sample={sample_name} {original_header} {read_number} {read_index}\n"
    return new_header


# Example usage
# process_fastq('/path/to/directory', sys.stdout)
# or to write to a file:
# with open('output.fastq', 'w') as output_file:
#     process_fastq('/path/to/directory', output_file)
