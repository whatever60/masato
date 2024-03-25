#!/usr/bin/env python3
import argparse
from io import StringIO
import sys
import subprocess
import os
import time
import tempfile

import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
from tqdm import tqdm, trange

from utils import print_command


def blast_online(
    input_path: str,
    database: str = "nt",
    batch_size: int = 100,
    email: str = "",
) -> list[str]:
    """
    Submit BLASTn jobs in batches to NCBI for sequences in a given FASTA file and save the results.

    Parameters:
    - input_path: Path to the input FASTA file containing the sequences to query.
    - database: The NCBI database to search against. Default is 'nt'.
    - output_results_path: Path where the combined BLAST results will be saved.
    - batch_size: The number of sequences to include in one query. Default is 100.
    """

    # Read sequences from the FASTA file
    sequences = list(SeqIO.parse(input_path, "fasta"))

    # Prepare batches using list slicing
    batches = [
        sequences[i : i + batch_size] for i in range(0, len(sequences), batch_size)
    ]

    results = []

    # Process each batch
    for batch in tqdm(batches, desc="Submitting BLAST queries"):
        # Convert batch to FASTA format string
        fasta_strings = "\n".join([record.format("fasta") for record in batch])

        result_handle = NCBIWWW.qblast(
            "blastn", database, fasta_strings, format_type="XML"
        )

        # Append result
        results.append(result_handle.read())
    return results


def parse_blast_results(result_dir: str) -> pd.DataFrame:
    results = []
    for xml_file in os.listdir(result_dir):
        if xml_file.endswith(".xml"):
            with open(os.path.join(result_dir, xml_file), "r") as f:
                results.append(f.read())
    return pd.concat([extract_info_from_xml(result, email) for result in results])

    # Write all results to the specified output file
    # with open(output_results_path, "w") as output_file:
    #     for result in results:
    #         output_file.write(result + "\n")

    # print(f"BLAST search completed and results saved to {output_results_path}")


def extract_info_from_xml(xml_string: str, email: str = "") -> pd.DataFrame:
    """
    Extract information from BLAST XML results and include taxonomy IDs.

    Parameters:
    - xml_string: BLAST results in XML format as a string.

    Returns:
    - A pandas DataFrame containing the results with taxonomy IDs.
    """
    records = []

    with StringIO(xml_string) as xml_handle:
        blast_records = NCBIXML.parse(xml_handle)
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    records.append(
                        dict(
                            qseqid=record.query,
                            qlen=record.query_length,
                            sseqid=alignment.hit_id,
                            pident=f"{hsp.identities / float(hsp.align_length) * 100:.2f}",
                            length=hsp.align_length,
                            qstart=hsp.query_start,
                            qend=hsp.query_end,
                            sstart=hsp.sbjct_start,
                            send=hsp.sbjct_end,
                            evalue=hsp.expect,
                            bitscore=hsp.bits,
                            slen=alignment.length,
                            saccession=alignment.accession,
                            sdefinition=alignment.hit_def,
                        )
                    )
    df = pd.DataFrame(records)
    acc2tax = fetch_taxonomy_ids(df["saccession"].unique(), email)
    df["staxids"] = df["saccession"].map(acc2tax)
    return df


def fetch_taxonomy_ids(
    accession_numbers: list[str], email: str = "", batch_size: int = 100
) -> dict:
    """
    Fetch taxonomy IDs for given accession numbers using NCBI's E-utilities in a single batch request.

    Parameters:
    - accession_numbers: A list of accession numbers.

    Returns:
    - A dictionary mapping accession numbers to taxonomy IDs.
    """
    Entrez.email = email
    taxonomy_ids = {}
    records = []
    # try:
    # Join accession numbers into a single string separated by commas
    for batch_idx, i in enumerate(trange(0, len(accession_numbers), batch_size)):
        successful = False
        while not successful:
            try:
                accession_str = ",".join(accession_numbers[i : i + batch_size])
                handle = Entrez.esummary(db="nucleotide", id=accession_str)
                records.extend(Entrez.read(handle))
                handle.close()
            except RuntimeError as e:
                print(f"An error occurred in batch {batch_idx}: {e}")
                print("Retrying in 5 seconds...")
                time.sleep(5)
            else:
                successful = True

    # Map accession numbers to taxonomy IDs
    for accession, record in zip(accession_numbers, records):
        # The 'Id' field contains the accession number; modify as needed based on the response structure
        tax_id = int(record["TaxId"])
        taxonomy_ids[accession] = tax_id
    # except Exception as e:
    #     print(f"An error occurred while fetching taxonomy IDs: {e}")

    return taxonomy_ids


def blast_local(input_path: str, database: str) -> pd.DataFrame:
    """performs blasts and returns the results"""
    # custom_blast_format = '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore slen staxids'
    custom_blast_format = "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore slen staxids"

    # db_command = ["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", "temp_db"] #command to 1econstruct database
    # sd = subprocess.Popen(db_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE); sd.communicate()
    # temp file
    temp_file = tempfile.NamedTemporaryFile()
    blast_command = [
        "megablast",
        "-query",
        input_path,
        "-db",
        database,
        "-num_threads",
        "16",
        "-evalue",
        "1e-30",
        "-out",
        temp_file.name,
        "-outfmt",
        custom_blast_format,
    ]
    print("Running blastn locally with command:")
    print_command(blast_command)
    sp = subprocess.Popen(
        blast_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )  # blast = sp.communicate()
    err = sp.stderr.read()
    if err:
        print(err)
        print("ERROR: with blast above")
        sys.exit()
    sp.communicate()
    blast_results = pd.read_table(
        temp_file.name,
        header=None,
        names=[
            "qseqid",
            "qlen",
            "sseqid",
            "pident",
            "length",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "slen",
            "staxids",
        ],
    )
    return blast_results


def merge_ranges(ranges: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping ranges and return the merged list."""
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    merged = []

    for start, end in sorted_ranges:
        if not merged or merged[-1][1] < start:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)

    return merged


def calculate_coverage(ranges: tuple[int, int], qlen: int) -> float:
    """Calculate the exact query coverage based on merged ranges and query length."""
    coverage = sum(end - start + 1 for start, end in ranges) / qlen
    return coverage * 100  # Convert to percentage


def add_query_coverage(df) -> pd.DataFrame:
    """Add query coverage to the DataFrame without modifying it in-place."""
    # Ensure not modifying the original DataFrame
    df_copy = df.copy()

    # Group by qseqid and aggregate qstart and qend into lists
    agg_df = (
        df_copy.groupby("qseqid")
        .agg({"qlen": "first", "qstart": list, "qend": list})
        .reset_index()
    )

    # Calculate merged ranges and coverage for each group
    agg_df["coverage"] = agg_df.apply(
        lambda x: calculate_coverage(
            merge_ranges(list(zip(x["qstart"], x["qend"]))), x["qlen"]
        ),
        axis=1,
    )

    # Merge the coverage back into the original DataFrame
    result_df = pd.merge(
        df_copy, agg_df[["qseqid", "coverage"]], on="qseqid", how="left"
    )

    return result_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input_",
        help="fasta file of sequences to be classified or a blast xml file",
    )
    parser.add_argument(
        "-o",
        "--blast_output",
        help="output file from blast",
    )
    parser.add_argument(
        "-d", "--database", help="fasta file of sequences to be classified"
    )
    parser.add_argument(
        "-b",
        "--batch_size",
        default=100,
        type=int,
        help="number of sequences to blast at once",
    )
    parser.add_argument("-e", "--email", help="email to use for NCBI E-utilities")

    args = parser.parse_args()

    input_ = args.input_
    database = args.database
    blast_output = args.blast_output
    batch_size = args.batch_size
    email = args.email
    if os.path.isdir(os.path.dirname(database)):
        df = blast_local(input_, database)
        df.to_csv(blast_output, index=False, sep="\t")
    elif os.path.isdir(input_):
        df = parse_blast_results(input_)
        df.to_csv(blast_output, index=False, sep="\t")
    else:
        # must be a fasta file
        results = blast_online(input_, database, batch_size, email)
        blast_result_dir = os.path.join(os.path.dirname(blast_output), "blast_results")
        for i, result in enumerate(results):
            with open(
                os.path.join(blast_result_dir, f"blast_result_{i}.xml"), "w"
            ) as f:
                print(result, file=f)
        df = parse_blast_results(blast_result_dir)
        df.to_csv(blast_output, index=False, sep="\t")
