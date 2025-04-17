import io
import glob
import gzip
import os
import re
import warnings
from typing import IO

import pandas as pd
import biom
from tqdm.auto import tqdm


PATTERN_ILLUMINA = re.compile(r"^(.+?)_S\d+_(?:L\d{3}_)?(R[12])_001.f(ast)?q(.gz)?$")
# - sample_R1.fq.gz or sample_R2.fq.gz (or fastq.gz, same for the following)
# - sample_1.fq.gz or sample_2.fq.gz
# - sample.R1.fq.gz or sample.R2.fq.gz
# - sample.1.fq.gz or sample.2.fq.gz
PATTERN_CUSTOM = re.compile(r"^(.+?)(?:(?:_|\.)(?:R)?[12])?.f(ast)?q(.gz)?$")


def find_paired_end_files(directory: str) -> list[tuple[str, str, str]]:
    # Dictionary to store file pairs
    file_pairs = {}

    # Find all .fastq.gz and .fq.gz files in the directory
    fastq_files = (
        glob.glob(os.path.join(directory, "*.fastq.gz"))
        + glob.glob(os.path.join(directory, "*.fq.gz"))
        + glob.glob(os.path.join(directory, "*.fastq"))
        + glob.glob(os.path.join(directory, "*.fq"))
    )

    # Process each file using pattern functions
    for file_ in fastq_files:
        for pattern in [PATTERN_ILLUMINA, PATTERN_CUSTOM]:
            match = pattern.search(os.path.basename(file_))
            if match:
                sample_name, read_type = match.groups()[:2]
                read_type = {"1": "R1", "2": "R2"}.get(read_type, read_type)
                if read_type not in ["R1", "R2"]:
                    warnings.warn(f"Invalid read type for file: {file_}")
                    continue

                pair_key = (sample_name, read_type)

                # Add the file to the dictionary and check for duplicates
                if pair_key in file_pairs:
                    warnings.warn(f"Duplicate file for {pair_key}: {file_}")
                else:
                    file_pairs[pair_key] = file_
                break
        else:
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


def smart_open(file_path: str, mode: str = "r") -> IO[str] | gzip.GzipFile:
    """Open a file as text or as a gzip file based on its magic number.

    Args:
        file_path (str): Path to the file to be opened.
        mode (str): Mode in which the file should be opened. Defaults to 'rt' (read text).

    Returns:
        IO[str] | gzip.GzipFile: A file object or a gzip file object.
    """

    if mode == "r":
        with open(file_path, "rb") as f:
            first_two_bytes = f.read(2)
        if first_two_bytes == b"\x1f\x8b":  # Magic number for gzip files
            return gzip.open(file_path, "rt")
        else:
            return open(file_path, "r")
    elif mode == "w":
        if file_path.endswith(".gz") or file_path.endswith(".gzip"):
            return gzip.open(file_path, "wt")
        else:
            return open(file_path, "w")
    else:
        raise ValueError(f"Mode must be 'r' or 'w' if specified, getting {mode}.")


def cat_fastq(
    directory: str,
    output_fp_r1,
    output_fp_r2=None,
    metadata: str | None = None,
    _remove_undet: bool = True,
    _have_sample_name: bool = False,
) -> None:
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
        samples_in_meta = None

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

    if _have_sample_name:
        rename_read = _rename_read_concat
    else:
        rename_read = _rename_read_illumina

    for idx, (r1_path, r2_path, sample_name) in enumerate(tqdm(matched_pairs)):
        if samples_in_meta is not None and sample_name not in samples_in_meta:
            continue
        if _remove_undet and sample_name == "Undetermined":
            continue
        # if idx < 4:
        #     continue
        # print(sample_name)
        with smart_open(r1_path) as r1_file, smart_open(r2_path) as r2_file:
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


def cat_fastq_se(
    directory: str,
    output_fp,
    metadata: str | None = None,
    _remove_undet: bool = True,
    _have_sample_name: bool = False,
    _r2: bool = False,
):
    """Similar as above but simply list all fastq/fq/fastq.gz/fq.gz files in the
    directory and concatenate them into a single file with new read names. Good for
    single-end reads.
    """
    fastq_files = (
        glob.glob(os.path.join(directory, "*.fastq.gz"))
        + glob.glob(os.path.join(directory, "*.fq.gz"))
        + glob.glob(os.path.join(directory, "*.fastq"))
        + glob.glob(os.path.join(directory, "*.fq"))
    )

    files = []
    for file_ in fastq_files:
        for pattern in [PATTERN_ILLUMINA, PATTERN_CUSTOM]:
            match = pattern.search(os.path.basename(file_))
            if match:
                sample_name, read_type = match.groups()[:2]
                read_type = {None: "R1", "1": "R1", "2": "R2"}.get(read_type, read_type)
                if read_type not in ["R1", "R2"]:
                    warnings.warn(f"Invalid read type {read_type} for file: {file_}")
                    continue
                if _r2:
                    if read_type == "R1":
                        break
                else:
                    if read_type == "R2":
                        break
                files.append((file_, sample_name, read_type))
                break

    if metadata is not None:
        samples_in_meta = pd.read_table(metadata, index_col=0).index.to_list()
    else:
        samples_in_meta = None

    if isinstance(output_fp, gzip.GzipFile) or isinstance(output_fp, io.BufferedWriter):

        def write_line(line):
            output_fp.write(line.encode())

    elif isinstance(output_fp, io.TextIOBase):

        def write_line(line):
            output_fp.write(line)

    else:
        raise ValueError("Output file pointer must be gzip or text file.")

    if _have_sample_name:
        rename_read = _rename_read_concat
    else:
        rename_read = _rename_read_illumina

    for file_, sample_name, read_type in tqdm(files):
        # match = PATTERN_ILLUMINA.search(os.path.basename(file_))
        # if match:
        # sample_name = match.group(1)
        # sample_name = os.path.basename(file_).split(".")[0].split("_")[0]
        if samples_in_meta is not None and sample_name not in samples_in_meta:
            continue
        if _remove_undet and sample_name == "Undetermined":
            continue
        with smart_open(file_) as f:
            for read_index, lines in enumerate(zip(*[f] * 4, strict=True), start=1):
                write_line(
                    rename_read(lines[0], sample_name, 1, read_index)
                    + "".join(lines[1:])
                )


def _rename_read_illumina(
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
    return f"@sample={sample_name} {read_number} {read_index} {original_header}\n"


def _rename_read_concat(
    header_line: str, sample_name: str, read_number: int, read_index: int
) -> str:
    header_line = header_line[1:]
    original_header, comment = header_line.split(maxsplit=1)
    original_sample = original_header.split("=", 1)[1]
    new_sample = f"{original_sample}_{sample_name}"
    return f"@sample={new_sample} {comment}"


def print_command(command: str | list[str]) -> None:
    """Print the given command in a readable format.
    Iterate through the command list and print each element on a new line if starts with
    a dash, otherwise print it on the same line.
    """
    print("Running the command:")
    if isinstance(command, str):
        print(command)
        return
    else:
        print(command[0], end="")
        for arg in command[1:]:
            if arg.startswith("-"):
                print(" \\\n    " + arg, end="")
            else:
                print(" " + arg, end="")
        print()


def read_table(
    table_path: str,
    index_col: str | int = 0,
    comment: str = None,
    dtype: str = "int",
    index_name: str = None,
) -> pd.DataFrame:
    """Read a table from a file and return it as a DataFrame.

    Args:
        table_path: Path to the table file. Could be a tsv, csv or biom.

    Returns:
        pd.DataFrame: DataFrame containing the table data.
    """
    if table_path.endswith(".biom"):
        df_biom = biom.load_table(table_path)
        df = df_biom.to_dataframe().astype(dtype)
        df.index.name = df_biom.table_id if index_name is None else index_name
        return df
    else:
        if table_path.endswith(".csv") or table_path.endswith(".csv.gz"):
            method = pd.read_csv
        elif table_path.endswith(".tsv") or table_path.endswith(".tsv.gz"):
            method = pd.read_table
        else:
            raise ValueError("Unsupported table file format.")
        return method(
            table_path, index_col=index_col, comment=comment, dtype={index_col: str}
        )


def write_table(table: pd.DataFrame, table_path: str) -> None:
    """Write a table to a file.

    Args:
        table: DataFrame to be written to a file.
        table_path: Path to the output file. Could be a tsv, csv or biom.
    """
    if table_path.endswith(".biom"):
        data = biom.Table(
            table.to_numpy(), table.index, table.columns, table_id=table.index.name
        )
        with biom.util.biom_open(table_path, "w") as f:
            data.to_hdf5(f, "whatever60")
    elif table_path.endswith(".csv") or table_path.endswith(".csv.gz"):
        table.to_csv(table_path)
    elif table_path.endswith(".tsv") or table_path.endswith(".tsv.gz"):
        table.to_csv(table_path, sep="\t")
    else:
        raise ValueError("Unsupported table file format.")
