import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"


def test_rdp(script_runner, test_data_dir):
    """
    run_rdp_classifier.py \
        -i $working_dir/$output_unoise3_dir/unoise3_zotu.fa \
        -d $rdp_db
    """
    # Prepare input and output directories
    input_file = test_data_dir / "output_unoise3_16s/unoise3_zotu.fa"
    db_file = "16srrna"

    # Run the script
    ret = script_runner.run(
        [
            "run_rdp_classifier.py",
            "-i",
            str(input_file),
            "-d",
            str(db_file),
        ]
    )

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"
    # Optionally, add more assertions to verify the content of the output file