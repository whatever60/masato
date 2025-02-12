#!/usr/bin/env python
import argparse
import subprocess
import os

import yaml

import masato


# RDP_CLASSIFIER_PATH: str = masato.__path__[0] + "/../../dist/rdp_classifier_2.14/dist/rdp_classifier"
# INFO_PATH = os.path.join(os.path.dirname(__file__), "setup/info.yaml")
# with open(INFO_PATH) as f:
#     info = yaml.safe_load(f)
#     rdp_classifier_path = info["rdp_classifier_path"]
#     if rdp_classifier_path is None:
#         RDP_CLASSIFIER_PATH = (
#             masato.__path__[0]
#             + "/../../dist/rdp_classifier_2.14/dist/rdp_classifier"
#         )
#     else:
#         RDP_CLASSIFIER_PATH = rdp_classifier_path
RDP_CLASSIFIER_PATH = (
    masato.__path__[0] + "/../../dist/rdp_classifier_2.14/dist/rdp_classifier"
)


def main():
    """
    Run these two commands with arguments for input zotu fasta and database. Output paths are inferred.

    rdp_classifier \
        $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
        -d $rdp_db \
        -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.tsv \
        -h $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.hier.tsv

    $src_dir/process_rrndb.py \
        -i $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb.tsv \
        -o $working_dir/$output_unoise3_dir/unoise3_zotu_rrndb_processed.tsv
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str)
    parser.add_argument("-d", "--database", type=str)

    args = parser.parse_args()
    input_path = args.input
    database_path = args.database
    base_name = os.path.splitext(input_path)[0]
    raw_output_path = f"{base_name}_rrndb_raw.tsv"
    raw_hier_path = f"{base_name}_rrndb_raw_hier.tsv"
    output_path = f"{base_name}_rrndb_processed.tsv"

    # Run the commands
    ret = subprocess.run(
        [
            RDP_CLASSIFIER_PATH,
            input_path,
            "-d",
            database_path,
            "-o",
            raw_output_path,
            "-h",
            raw_hier_path,
        ]
    )
    assert (
        ret.returncode == 0
    ), f"rdp_classifier failed with return code {ret.returncode}"

    ret = subprocess.run(["process_rrndb.py", "-i", raw_output_path, "-o", output_path])
    assert (
        ret.returncode == 0
    ), f"process_rrndb.py failed with return code {ret.returncode}"
