#!/usr/bin/python3
import argparse

import pandas as pd
from tqdm.auto import tqdm

# Create the parser
parser = argparse.ArgumentParser(description="Process names and nodes files.")

# Add the arguments
parser.add_argument("names_file", type=str, help="The names file")
parser.add_argument("nodes_file", type=str, help="The nodes file")

# Parse the arguments
args = parser.parse_args()

names_file = args.names_file
nodes_file = args.nodes_file

# master_dict[id] = [rank,name_text,parent_tax_id]

print("Parsing names.dmp file for tax_id and scientific_name...")
# Read names file into DataFrame
# tax_id: name_text
names_dictionary = (
    pd.read_csv(
        names_file,
        sep="|",
        header=None,
        usecols=[0, 1, 2, 3],
        names=["tax_id", "name_text", "unique_name", "name_class"],
        dtype=str,
    )
    .map(lambda x: x.strip() if isinstance(x, str) else x)
    .query('name_class == "scientific name"')
    .drop_duplicates("tax_id")
    .set_index("tax_id")["name_text"]
    .to_dict()
)
print(f"A total of {len(names_dictionary)} unique entries in the names database")

print("Parsing nodes.dmp for tax_id, parent_tax_id and taxonomy rank...")
# Read nodes file into DataFrame
nodes_dictionary = (
    pd.read_csv(
        nodes_file,
        sep="|",
        header=None,
        usecols=[0, 1, 2],
        names=["tax_id", "parent_tax_id", "rank"],
        dtype=str,
    )
    .map(lambda x: x.strip() if isinstance(x, str) else x)
    .set_index("tax_id")[["rank", "parent_tax_id"]]
    .transpose()
    .to_dict(orient="list")
)
# tax_id: rank, parent_tax_id
print(f"A total of {len(nodes_dictionary)} unique entries in the names database")

import pdb; pdb.set_trace()
def expand_taxonomy(tax_id, names_dictionary, nodes_dictionary):
    tax_levels = [
        "superkingdom",
        "kingdom",
        "phylum",
        "subphylum",
        "superclass",
        "class",
        "subclass",
        "superorder",
        "order",
        "superfamily",
        "family",
        "subfamily",
        "genus",
        "species",
    ]
    tax_level2idx = {tax_levels[i]: i for i in range(len(tax_levels))}
    default_tax = ["unknown_" + x for x in tax_levels]
    growing_taxonomy = [names_dictionary[tax_id]]
    growing_taxonomy_levels = [nodes_dictionary[tax_id][0]]
    master_key = tax_id

    while True:
        new_tax_id = nodes_dictionary[tax_id][1]  # parent_tax_id
        # new_name = names_dictionary[new_tax_id]
        # new_rank = nodes_dictionary[new_tax_id][0]
        new_name, new_rank = (
            names_dictionary[new_tax_id],
            nodes_dictionary[new_tax_id][0],
        )
        if new_name == "root":
            break
        growing_taxonomy.append(new_name)
        growing_taxonomy_levels.append(new_rank)
        tax_id = new_tax_id
    reversed_taxonomy = list(reversed(growing_taxonomy))
    reversed_taxonomy_levels = list(reversed(growing_taxonomy_levels))

    # final_tax = ['unknown_superkingdom', 'unknown_kingdom', 'unknown_phylum', 'unknown_subphylum', 'unknown_superclass', 'unknown_class', 'unknown_subclass', 'unknown_superorder', 'unknown_order', 'unknown_superfamily', 'unknown_family', 'unknown_subfamily', 'unknown_family', 'unknown_genus', 'unknown_species']

    # for t in range(len(tax_levels)):
    #     for i in range(len(reversed_taxonomy_levels)):
    #         if reversed_taxonomy_levels[i] == tax_levels[t]:
    #             final_tax[t] = reversed_taxonomy[i]
    for l, t in zip(reversed_taxonomy_levels, reversed_taxonomy):
        if l in tax_level2idx:
            default_tax[tax_level2idx[l]] = t

    # format_tax = []
    # for i in range(len(final_tax)):
    #     format_tax.append('D_'+str(i)+'__'+final_tax[i])
    # combined_taxonomy = ';'.join(format_tax)

    # print (master_key+'\t'+combined_taxonomy)
    return master_key, ";".join(default_tax)


output_file = open("expanded_ncbi_taxonomy.tsv", "w")

for tax_id in tqdm(names_dictionary):  # .sort(key=int):
    t, expanded_taxonomy = expand_taxonomy(tax_id, names_dictionary, nodes_dictionary)
    output_file.writelines(t + "\t" + expanded_taxonomy + "\n")
