import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"


@pytest.mark.parametrize("amplicon_type, prefix", [("16s", "16S-U")])
def test_unoise3(script_runner, test_data_dir, amplicon_type, prefix):
    """
    Test for running the unoise3 command.
    """
    # Prepare input and output directories
    output_unoise3_dir = f"output_unoise3_{amplicon_type}"
    input_file = test_data_dir / output_unoise3_dir / "merged.fq.gz"
    output_file = test_data_dir / output_unoise3_dir / "unoise3_zotu.fa"

    # Run the script
    ret = script_runner.run(
        [
            "usearch_workflow.py",
            "unoise3",
            "-i",
            str(input_file),
            "-o",
            str(output_file),
            "-l",
            prefix,
        ]
    )

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"
    assert output_file.exists(), f"Output file not created: {output_file}"
    # Optionally, add more assertions to verify the content of the output file


@pytest.mark.parametrize(
    "amplicon_type, unknown_name", [("16s", "16S-U_UNKNOWN")]
)
def test_search_global(
    script_runner, test_data_dir, amplicon_type, unknown_name
):
    """
    Test for running the search_global command.
    """
    # Prepare input and output directories
    output_unoise3_dir = f"output_unoise3_{amplicon_type}"
    input_file = test_data_dir / output_unoise3_dir / "merged.fq.gz"
    db_file = test_data_dir / output_unoise3_dir / "unoise3_zotu.fa"
    output_tsv = test_data_dir / output_unoise3_dir / "unoise3_zotu.biom"

    # Run the script
    ret = script_runner.run(
        [
            "usearch_workflow.py",
            "search_global",
            "-i",
            str(input_file),
            "-d",
            str(db_file),
            "-o",
            str(output_tsv),
            "--unknown_name",
            unknown_name,
        ]
    )

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"
    assert output_tsv.exists(), f"Output file not created: {output_tsv}"
    # Optionally, add more assertions to verify the content of the output file


@pytest.mark.parametrize("amplicon_type, prefix", [("16s", "16S-U")])
def test_usearch_workflow(script_runner, test_data_dir, amplicon_type, prefix):
    """
    Test script for running two `usearch_workflow.py` commands:
    1. `workflow_per_sample`
    2. `aggregate_samples`
    """

    # Define test variables
    output_unoise3_dir = f"output_unoise3_{amplicon_type}"

    # Prepare input and output directories
    input_file = test_data_dir / output_unoise3_dir / "merged_isolate.fq.gz"
    json_output_file = test_data_dir / output_unoise3_dir / "unoise3_zotu_isolate.json"
    fasta_output_file = test_data_dir / output_unoise3_dir / "unoise3_zotu_isolate.fa"
    biom_output_file = test_data_dir / output_unoise3_dir / "unoise3_zotu_isolate.biom"

    # First command: workflow_per_sample
    cmd_workflow_per_sample = [
        "usearch_workflow.py",
        "workflow_per_sample",
        "-i",
        str(input_file),
        "-o",
        str(json_output_file),
        # "-l",
        # "whatever",
        "--search",
    ]

    # Run the first command
    ret_workflow_per_sample = script_runner.run(cmd_workflow_per_sample)

    # Assertions for the first command
    assert (
        ret_workflow_per_sample.success
    ), f"Script failed with return code {ret_workflow_per_sample.returncode}"
    assert json_output_file.exists(), f"Output file not created: {json_output_file}"

    # Second command: aggregate_samples
    cmd_aggregate_samples = [
        "usearch_workflow.py",
        "aggregate_samples",
        "-i",
        str(json_output_file),
        "-of",
        str(fasta_output_file),
        "-oc",
        str(biom_output_file),
        "-p",
        prefix,
    ]

    # Run the second command
    ret_aggregate_samples = script_runner.run(*cmd_aggregate_samples)

    # Assertions for the second command
    assert (
        ret_aggregate_samples.success
    ), f"Script failed with return code {ret_aggregate_samples.returncode}"
    assert fasta_output_file.exists(), f"Output file not created: {fasta_output_file}"
    assert biom_output_file.exists(), f"Output file not created: {biom_output_file}"
