import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"


def test_get_tree(script_runner, test_data_dir):
    """
    $src_dir/get_tree.py -i $working_dir/$output_unoise3_dir/unoise3_zotu.fa
    """
    # Prepare input and output directories
    input_file = test_data_dir / "output_unoise3_16s/unoise3_zotu.fa"
    output_file = test_data_dir / "output_unoise3_16s/unoise3_zotu.newick"

    # Run the script
    ret = script_runner.run(["get_tree.py", "-i", str(input_file)])

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"
    assert output_file.exists(), f"Output file not created: {output_file}"
    # Optionally, add more assertions to verify the content of the output file
