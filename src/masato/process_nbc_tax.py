"""Parse the output of usearch -nbc_tax to extract taxonomic classification.

Two tsv files will be generated, one representing complete information in the input 
file and another only extract annotation at genus and family level. For the latter, 
any OTUs without assignment at genus or family level are discarded.
"""

import argparse
import re

import pandas as pd

print(
    "WARNING: This script is deprecated since it is for processing the output of "
    "`usearch -nbc_tax` command (which reimplements RPD classifier) and `usearch` is "
    "no longer used in this project. This script is here just for bookkeeping purpose. "
    "Use the original RDP classifier JAVA implementation instead."
)


def parse_taxonomy_line(row):
    # Define regex patterns
    prob_pattern = re.compile(r"([a-z]):([^,:]+)\((\d+\.\d+)\)")
    assign_pattern = re.compile(r"([a-z]):([^,:]+)")

    prob_matches = prob_pattern.findall(row["prob"])
    assign_matches = (
        assign_pattern.findall(row["assign"]) if pd.notna(row["assign"]) else []
    )

    # Create an empty result dictionary
    result = {
        "otu": row["otu"],
        "domain": None,
        "domain_p": None,
        "domain_conf": None,
        "phylum": None,
        "phylum_p": None,
        "phylum_conf": None,
        "class": None,
        "class_p": None,
        "class_conf": None,
        "order": None,
        "order_p": None,
        "order_conf": None,
        "family": None,
        "family_p": None,
        "family_conf": None,
        "genus": None,
        "genus_p": None,
        "genus_conf": None,
        "species": None,
        "species_p": None,
        "species_conf": None,
    }

    # Mapping for the taxonomic ranks
    tax_mapping = {
        "d": "domain",
        "p": "phylum",
        "c": "class",
        "o": "order",
        "f": "family",
        "g": "genus",
        "s": "species",
    }

    # Populate the result dictionary with classification matches
    for rank, tax, prob in prob_matches:
        if rank in tax_mapping:
            result[tax_mapping[rank]] = tax
            result[tax_mapping[rank] + "_p"] = float(prob)
    # Populate the result dictionary with confidence based on assign
    for rank, tax in assign_matches:
        if rank in tax_mapping and result[tax_mapping[rank]] == tax:
            result[tax_mapping[rank] + "_conf"] = tax

    return pd.Series(result)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str)
    parser.add_argument("output", type=str)
    parser.add_argument("output_simple", type=str)

    args = parser.parse_args()
    input_path = args.input
    output_path = args.output
    output_simple_path = args.output_simple

    # Load the data from the tsv
    df = pd.read_csv(
        input_path,
        sep="\t",
        names=["otu", "prob", "strand", "assign"],
    )

    # Process each row
    df_res = df.apply(parse_taxonomy_line, axis=1)
    df_res.to_csv(
        output_path,
        sep="\t",
        index=False,
    )

    df_res = df_res[["otu", "genus_conf", "family_conf"]].dropna()
    df_res.columns = ["otu", "genus", "family"]
    df_res.to_csv(
        output_simple_path,
        sep="\t",
        index=False,
    )
