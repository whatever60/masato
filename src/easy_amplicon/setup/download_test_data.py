#!/usr/bin/env python
"""
Simply run wget -qO- http://easy-amplicon-camii-test-data.s3.amazonaws.com/data.tar.gz | tar xz -C <user_specified_dir>
"""
import argparse
import subprocess


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--date", type=str)
    parser.add_argument("-o", "--output_dir", type=str)
    args = parser.parse_args()
    date = args.date
    output_dir = args.output_dir
    res = subprocess.run(
        f"wget -qO- http://easy-amplicon-camii-test-data.s3.amazonaws.com/data_{date}.tar.gz | tar xz -C {output_dir}",
        shell=True,
        check=True,
    )
    assert res.returncode == 0, f"Download failed with return code {res.returncode}"
