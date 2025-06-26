from calendar import c
import io
import glob
import gzip
import sys
import os
import re
import warnings
from typing import IO, BinaryIO, TextIO, Callable

import pandas as pd
from tqdm.auto import tqdm

from .io import read_table, write_table  # used by other modules.


PATTERN_ILLUMINA = re.compile(r"^(.+?)_S\d+_(?:L\d{3}_)?(R[12])_001.f(ast)?q(.gz)?$")
# - sample_R1.fq.gz or sample_R2.fq.gz (or fastq.gz, same for the following)
# - sample_1.fq.gz or sample_2.fq.gz
# - sample.R1.fq.gz or sample.R2.fq.gz
# - sample.1.fq.gz or sample.2.fq.gz
PATTERN_CUSTOM = re.compile(r"^(.+?)(?:(?:_|\.)(R?[12]))?.f(?:ast)?q(?:.gz)?$")


def find_paired_end_files(fastq_input: str | list[str]) -> list[tuple[str, str, str]]:
    """Identify matched paired-end FASTQ files based on filename patterns.

    Supports input as a directory path, a glob pattern (e.g., 'data/*.fastq.gz'),
    or a list of FASTQ file paths. Matches standard Illumina naming conventions
    and custom formats like sample_R1.fq.gz, sample.2.fastq, etc.

    Args:
        fastq_input (str | list[str]): A path to a directory, a glob pattern string,
            or a list of FASTQ file paths.

    Returns:
        list[tuple[str, str, str]]: A list of (R1_path, R2_path, sample_name) tuples,
            sorted by sample name.

    Raises:
        ValueError: If no FASTQ files are found or if no paired-end files can be matched.
    """
    if isinstance(fastq_input, str):
        if os.path.isdir(fastq_input):
            fastq_files = (
                glob.glob(os.path.join(fastq_input, "*.fastq.gz"))
                + glob.glob(os.path.join(fastq_input, "*.fq.gz"))
                + glob.glob(os.path.join(fastq_input, "*.fastq"))
                + glob.glob(os.path.join(fastq_input, "*.fq"))
            )
        else:
            fastq_files = glob.glob(fastq_input)
    else:
        fastq_files = fastq_input

    if not fastq_files:
        raise ValueError("No FASTQ files found in the specified input.")

    file_pairs: dict[tuple[str, str], str] = {}
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

    if not matched_pairs:
        raise ValueError("No paired-end FASTQ file pairs were found.")

    # Report the number of matched pairs
    print(f"Number of matched pairs: {len(matched_pairs)}")

    # Order by sample name
    matched_pairs.sort(key=lambda x: x[2])
    return matched_pairs


def smart_open(file_path: str, mode: str = "r") -> io.TextIOWrapper | gzip.GzipFile:
    """Open a file as text or as a gzip file based on its magic number.

    Args:
        file_path (str): Path to the file to be opened.
        mode (str): Mode in which the file should be opened. Defaults to 'rt' (read text).

    Returns:
        IO[str] | IO[byte]: A file object or a gzip file object.
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


OutputHandle = str | None | TextIO | BinaryIO


def resolve_output(fp: OutputHandle) -> IO[str] | IO[bytes]:
    if fp is None or fp == "-":
        return sys.stdout
    elif isinstance(fp, (io.TextIOBase, io.StringIO, io.BufferedIOBase, io.BytesIO)):
        return fp
    elif isinstance(fp, str):
        return smart_open(fp, mode="w")
    else:
        raise ValueError(f"Unsupported output type: {type(fp)}")


def _make_writer(handle: IO[str] | IO[bytes]) -> Callable[[str], None]:
    """
    Create a string-writing function based on the stream type.

    This function returns a writer that takes a string and writes it to the given handle.
    If the handle is a text stream (`io.TextIOBase`), the string is written directly.
    If the handle is a binary stream (`io.BufferedIOBase`), the string is encoded to UTF-8
    bytes before writing.

    Common use cases:
        - `sys.stdout`: Usually `TextIOBase`; safe to write strings directly.
        - `sys.stdout.buffer`: `BufferedIOBase`; requires encoding.
        - `open("file.txt", "w")`: Returns a `TextIOWrapper` (text mode).
        - `open("file.txt", "wb")`: Returns a `BufferedWriter` (binary mode).
        - `subprocess.Popen(..., stdin=PIPE, text=True)`: `proc.stdin` is `TextIOBase`.
        - `subprocess.Popen(..., stdin=PIPE, text=False)`: `proc.stdin` is `BufferedIOBase`.
        - `io.StringIO()`: `TextIOBase`, accepts strings.
        - `io.BytesIO()`: `BufferedIOBase`, accepts bytes.

    Args:
        handle: A file-like object for writing, either in text mode (`IO[str]`) or binary mode (`IO[bytes]`).

    Returns:
        A function that takes a string and writes it to the handle appropriately.

    Raises:
        ValueError: If the stream type is not recognized as either text or binary.
    """

    if isinstance(handle, io.TextIOBase):

        def func(s: str) -> None:
            handle.write(s)
    elif isinstance(handle, io.BufferedIOBase):

        def func(s: str) -> None:
            handle.write(s.encode())
    else:
        raise ValueError("Unsupported stream type")
    return func


def cat_fastq(
    directory: str | list[str],
    output_fp_r1: OutputHandle,
    output_fp_r2: OutputHandle = None,
    metadata: str | None = None,
    _remove_undet: bool = True,
    _have_sample_name: bool = False,
) -> None:
    """
    Concatenate paired-end FASTQ files from a directory into one or two output streams,
    renaming reads along the way. Supports output to files, subprocess pipes, stdout,
    or interleaved output depending on destination identity.

    Args:
        directory: Directory containing FASTQ files to process.
        output_fp_r1: Destination for R1 reads. Can be a file path, a writable file-like
            object (text or binary), subprocess pipe (`proc.stdin`), or "-" / None for stdout.
        output_fp_r2: Destination for R2 reads. Same options as `output_fp_r1`. If None,
            it defaults to stdout. If both outputs are the same object, interleaved output
            is written.
        metadata: Optional path to a sample metadata file. If provided, only samples listed
            in this file (1st column as index) are processed.
        _remove_undet: Whether to skip samples named "Undetermined".
        _have_sample_name: If True, include sample name in renamed read ID. Otherwise use
            standard Illumina-style renaming.

    Raises:
        ValueError: If output types are unsupported.
    """
    if metadata is not None:
        samples_in_meta = pd.read_table(metadata, index_col=0).index.to_list()
    else:
        samples_in_meta = None

    out_r1 = resolve_output(output_fp_r1)
    out_r2 = resolve_output(output_fp_r2)

    write_r1 = _make_writer(out_r1)
    write_r2 = _make_writer(out_r2)

    rename_read = _rename_read_concat if _have_sample_name else _rename_read_illumina

    for idx, (r1_path, r2_path, sample_name) in enumerate(
        tqdm(find_paired_end_files(directory))
    ):
        if samples_in_meta is not None and sample_name not in samples_in_meta:
            continue
        if _remove_undet and sample_name == "Undetermined":
            continue

        with smart_open(r1_path) as r1_file, smart_open(r2_path) as r2_file:
            paired_read_iter = zip(
                zip(*[r1_file] * 4, strict=True),
                zip(*[r2_file] * 4, strict=True),
                strict=True,
            )
            for read_index, (r1_lines, r2_lines) in enumerate(
                paired_read_iter, start=1
            ):
                r1_block = rename_read(
                    r1_lines[0], sample_name, 1, read_index
                ) + "".join(r1_lines[1:])
                r2_block = rename_read(
                    r2_lines[0], sample_name, 2, read_index
                ) + "".join(r2_lines[1:])
                write_r1(r1_block)
                write_r2(r2_block)


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

    prog_bar = tqdm(files)
    for file_, sample_name, read_type in prog_bar:
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
    original_header = header_line[1:].strip()
    return f"@{original_header} {read_index} ;sample={sample_name}\n"


def _rename_read_concat(
    header_line: str, sample_name: str, read_number: int, read_index: int
) -> str:
    original_header = header_line[1:].strip()
    assert original_header.split()[-1].startswith(";sample="), (
        "Header line must end with ';sample='"
    )
    return f"@{original_header}_{sample_name}\n"


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


# def read_table(
#     table_path: str,
#     index_col: str | int = 0,
#     comment: str = None,
#     dtype: str = "int",
#     index_name: str = None,
# ) -> pd.DataFrame:
#     """Read a table from a file and return it as a DataFrame.

#     Args:
#         table_path: Path to the table file. Could be a tsv, csv or biom.

#     Returns:
#         pd.DataFrame: DataFrame containing the table data.
#     """
#     if table_path.endswith(".biom"):
#         df_biom = biom.load_table(table_path)
#         df = df_biom.to_dataframe().astype(dtype)
#         df.index.name = df_biom.table_id if index_name is None else index_name
#         return df
#     else:
#         if table_path.endswith(".csv") or table_path.endswith(".csv.gz"):
#             method = pd.read_csv
#         elif table_path.endswith(".tsv") or table_path.endswith(".tsv.gz"):
#             method = pd.read_table
#         else:
#             raise ValueError("Unsupported table file format.")
#         return method(
#             table_path, index_col=index_col, comment=comment, dtype={index_col: str}
#         )


# def write_table(table: pd.DataFrame, table_path: str) -> None:
#     """Write a table to a file.

#     Args:
#         table: DataFrame to be written to a file.
#         table_path: Path to the output file. Could be a tsv, csv or biom.
#     """
#     if table_path.endswith(".biom"):
#         data = biom.Table(
#             table.to_numpy(), table.index, table.columns, table_id=table.index.name
#         )
#         with biom.util.biom_open(table_path, "w") as f:
#             data.to_hdf5(f, "whatever60")
#     elif table_path.endswith(".csv") or table_path.endswith(".csv.gz"):
#         table.to_csv(table_path)
#     elif table_path.endswith(".tsv") or table_path.endswith(".tsv.gz"):
#         table.to_csv(table_path, sep="\t")
#     else:
#         raise ValueError("Unsupported table file format.")
