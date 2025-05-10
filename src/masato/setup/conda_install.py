#!/usr/bin/env python
"""
Simply run conda install --file conda_requirements.txt -c defaults -c bioconda -c conda-forge,
where the requirement file is in the root of the package. (This file in /src/masato/setup/conda_install.py)
"""

import argparse
import subprocess
import os

from ..utils import print_command


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--flavor", type=str, default="conda", choices=["conda", "mamba"]
    )
    args = parser.parse_args()
    flavor = args.flavor
    requirement_file = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "data", "conda_requirements.txt")
    )
    args = [
        flavor,
        "install",
        "--file",
        requirement_file,
        "-c",
        "defaults",
        "-c",
        "bioconda",
        "-c",
        "conda-forge",
        "-y",
    ]
    print_command(args)
    ret = subprocess.run(args)
    assert ret.returncode == 0, (
        f"Conda install failed with return code {ret.returncode}"
    )
