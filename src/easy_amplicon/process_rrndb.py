#!/usr/bin/env python3

import argparse

import pandas as pd


def parse_taxonomy_line(row):
    # Prepare the taxonomy and probability entries
    start_index = 5
    entries = iter(row[start_index:])
    # prob_entries = list(row[start_index::3])
    # tax_entries = list(row[start_index + 1 :: 3])
    # rank_entries = list(row[start_index + 2 :: 3])

    # Mapping for the taxonomic ranks
    tax_mapping = {
        "domain": "domain",
        "phylum": "phylum",
        "class": "class",
        "order": "order",
        "family": "family",
        "genus": "genus",
        "species": "species",
    }

    # Create an empty result dictionary
    result = {
        "otu": row[0],
        "domain": "unknown",
        "domain_p": -1,
        "phylum": "unknown",
        "phylum_p": -1,
        "class": "unknown",
        "class_p": -1,
        "order": "unknown",
        "order_p": -1,
        "family": "unknown",
        "family_p": -1,
        "genus": "unknown",
        "genus_p": -1,
        "species": "unknown",
        "species_p": -1,
    }

    # Populate the result dictionary with the provided data
    taxon_names = set()
    for tax, rank, prob in zip(entries, entries, entries):
        if rank in tax_mapping:
            if tax != "unknwon" and tax in taxon_names:
                tax += rank[0]
            tax = tax.replace(" ", "_")
            taxon_names.add(tax)
            result[tax_mapping[rank]] = tax
            result[tax_mapping[rank] + "_p"] = float(prob)

    return pd.Series(result)


def main():
# if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str)
    parser.add_argument("-o", "--output", type=str)

    args = parser.parse_args()
    input_path = args.input
    output_path = args.output

    # Load the data from the tsv
    df = pd.read_table(input_path, header=None, names=list(range(30)))

    # Process each row
    df_res = df.apply(parse_taxonomy_line, axis=1)
    df_res.to_csv(output_path, sep="\t", index=False)
