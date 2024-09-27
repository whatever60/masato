import subprocess
from pathlib import Path

import pytest
from pytest_console_scripts import ScriptRunner


@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"


@pytest.fixture
def barcodes_files():
    return {
        "barcodes_fwd": "isolate/barcodes_fwd.fasta",
        "barcodes_rev": "isolate/barcodes_rev.fasta",
        "patterns": "isolate/patterns.txt",
    }


@pytest.mark.parametrize(
    "mode, input_folder, output_file_name, use_barcodes",
    [
        ("isolate_150", "fastq_isolate_16s", "merged_isolate.fq.gz", True),
        ("r1", "fastq_bulk_16s", "merged_bulk.fq.gz", False),
        ("simple", "fastq_bulk_16s", "merged.fq.gz", False),  # Added the simple mode
    ],
)
def test_trim(
    script_runner: ScriptRunner,
    test_data_dir,
    barcodes_files,
    mode,
    input_folder,
    output_file_name,
    use_barcodes,
):
    """
    Generalized test for `trim.py` with three modes: `isolate_150`, `r1`, and `simple`.
    """
    input_dir = test_data_dir / input_folder
    output_file = test_data_dir / "output_unoise3_16s" / output_file_name
    output_file.parent.mkdir(parents=True, exist_ok=True)

    # If the mode is `isolate_150`, we will pass barcodes
    extra_args = []
    if use_barcodes:
        barcodes_fwd = barcodes_files["barcodes_fwd"]
        barcodes_rev = barcodes_files["barcodes_rev"]
        patterns = barcodes_files["patterns"]
        extra_args = ["-fb", barcodes_fwd, "-rb", barcodes_rev, "-pt", patterns]

    # Run the script
    ret = script_runner.run(
        [
            "trim.py",
            "--mode",
            mode,
            "-i",
            str(input_dir),
            "-o",
            str(output_file),
            "-p",
            "16s",
            "-k",
            "110",
            "-l",
            "110",
        ]
        + extra_args
    )

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"
    assert output_file.exists(), f"Output file not created: {output_file}"
    # Optionally, add more assertions to verify the content of the output file


# def test_merge(test_data_dir):
#     """
#     Merge `merged_bulk.fq.gz` and `merged_isolate.fq.gz` using `zcat` and `gzip`.
#     """
#     bulk_file = test_data_dir / "output_unoise3_16s" / "merged_bulk.fq.gz"
#     isolate_file = test_data_dir / "output_unoise3_16s" / "merged_isolate.fq.gz"
#     output_file = test_data_dir / "output_unoise3_16s" / "merged.fq.gz"

#     command = f"zcat {bulk_file} {isolate_file} | gzip > {output_file}"

#     # Run the shell command
#     proc = subprocess.Popen(
#         command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
#     )
#     stdout, stderr = proc.communicate()

#     # Assertions
#     assert proc.returncode == 0, f"Script failed with return code {proc.returncode}. Error: {stderr.decode()}"
#     assert output_file.exists(), f"Output file not created: {output_file}"
