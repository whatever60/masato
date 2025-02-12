#!/usr/bin/env python
"""
Simply run wget -qO- http://easy-amplicon-camii-test-data.s3.amazonaws.com/dist.tar.gz | tar xz -C masato.__path__[0] + "/../.."
"""

import subprocess

import masato


def main():
    output_dir = masato.__path__[0] + "/../.."

    res = subprocess.run(
        f"wget -qO- http://easy-amplicon-camii-test-data.s3.amazonaws.com/dist.tar.gz | tar xz -C {output_dir}", shell=True, check=True
    )
    assert res.returncode == 0, f"Download failed with return code {res.returncode}"

    # write to the info yaml file {output_dir}/rdp_classifier_2.14 to the rdp_classifier_path key
    # with open(INFO_PATH) as f:
    #     info = yaml.safe_load(f)
    # info["rdp_classifier_path"] = f"{output_dir}/rdp_classifier_2.14"
    # with open(INFO_PATH, "w") as f:
    #     yaml.dump(info, f)
