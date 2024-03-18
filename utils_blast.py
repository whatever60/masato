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
from tqdm import tqdm

from utils import print_command


def blast_online(
    input_fasta_path: str,
    database: str = "nt",
    batch_size: int = 100,
    email: str = "",
) -> pd.DataFrame:
    """
    Submit BLASTn jobs in batches to NCBI for sequences in a given FASTA file and save the results.

    Parameters:
    - input_fasta_path: Path to the input FASTA file containing the sequences to query.
    - database: The NCBI database to search against. Default is 'nt'.
    - output_results_path: Path where the combined BLAST results will be saved.
    - batch_size: The number of sequences to include in one query. Default is 100.
    """

    # Read sequences from the FASTA file
    sequences = list(SeqIO.parse(input_fasta_path, "fasta"))

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
    import pdb; pdb.set_trace()
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


def fetch_taxonomy_ids(accession_numbers: list[str], email: str = "") -> dict:
    """
    Fetch taxonomy IDs for given accession numbers using NCBI's E-utilities in a single batch request.

    Parameters:
    - accession_numbers: A list of accession numbers.

    Returns:
    - A dictionary mapping accession numbers to taxonomy IDs.
    """
    Entrez.email = email
    taxonomy_ids = {}
    try:
        # Join accession numbers into a single string separated by commas
        accession_str = ",".join(accession_numbers)
        handle = Entrez.esummary(db="nucleotide", id=accession_str)
        records = Entrez.read(handle)
        handle.close()

        # Map accession numbers to taxonomy IDs
        for record in records:
            # The 'Id' field contains the accession number; modify as needed based on the response structure
            accession = record["AccessionVersion"].split(".")[0]
            tax_id = record["TaxId"]
            taxonomy_ids[accession] = tax_id
    except Exception as e:
        print(f"An error occurred while fetching taxonomy IDs: {e}")

    return taxonomy_ids


def blast_local(input_fasta_path: str, database: str) -> pd.DataFrame:
    """performs blasts and returns the results"""
    # custom_blast_format = '6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore slen staxids'
    custom_blast_format = "6 qseqid qlen sseqid pident length qstart qend sstart send evalue bitscore slen staxids"

    # db_command = ["makeblastdb", "-in", subject, "-dbtype", "nucl", "-out", "temp_db"] #command to 1econstruct database
    # sd = subprocess.Popen(db_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE); sd.communicate()
    # temp file
    temp_file = tempfile.NamedTemporaryFile()
    blast_command = [
        "blastn",
        "-query",
        input_fasta_path,
        "-db",
        database,
        "-num_threads",
        "16",
        "-evalue",
        "1e-5",
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
        "-i", "--input_fasta", help="fasta file of sequences to be classified"
    )
    # parser.add_argument(
    #     "-n",
    #     "--ncbi_nt",
    #     help="use the NCBI NT database for taxonomy",
    #     action="store_true",
    # )
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

    # if args.ncbi_nt:
    #     ###INFO FOR TAXONOMY DATABASES
    #     taxonomy_categories = [
    #         "superkingdom",
    #         "kingdom",
    #         "phylum",
    #         "subphylum",
    #         "superclass",
    #         "class",
    #         "subclass",
    #         "superorder",
    #         "order",
    #         "superfamily",
    #         "family",
    #         "subfamily",
    #         "genus",
    #         "species",
    #     ]
    #     # qiime_all_level_12_categories = [superkingdom, subkingdom, sub_subkingdom, kingdom, tmp1, tmp2, phylum, subphylum, class, order, family, genus, species]
    #     seven_levels = [0, 1, 2, 5, 8, 10, 12, 13]
    #     ### OTHER OPTIONS BASED ON DATABASE OF CHOICE
    #     uncultured_cutoff_level = len(taxonomy_categories)
    #     species_level = 13
    #     family_level = 10
    #     phylum_level = 2

    #     taxonomy_to_grab = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]

    # else:
    #     seven_levels = [0, 1, 2, 5, 8, 10, 12, 13]

    # sequence_dict = SeqIO.to_dict(SeqIO.parse(args.input_fasta, "fasta"))
    # taxonomy_dictionary = {}  # tax_code:[taxonomy]
    # sequence_taxonomy_dict = {}  # header:[best_taxonomy]
    # percent_id_dict = {}
    input_fasta = args.input_fasta
    database = args.database
    blast_output = args.blast_output
    batch_size = args.batch_size
    email = args.email
    if os.path.isdir(os.path.dirname(database)):
        df = blast_local(input_fasta, database)
    else:
        df = blast_online(input_fasta, database, batch_size, email)
    df.to_csv(blast_output, index=False, sep="\t")
