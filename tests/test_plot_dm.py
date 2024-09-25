import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"


@pytest.mark.parametrize("amplicon_type", ["16s"])
@pytest.mark.parametrize("level", ["genus"])
@pytest.mark.parametrize(
    "distance_method, log10, annotated, suffix",
    [
        ("euclid", True, False, ""),  # Euclidean with log10, no annotation
        ("braycurtis", False, False, ""),  # Bray-Curtis without log10, no annotation
        ("euclid", True, True, "annot_"),  # Euclidean with log10 and annotation
        ("braycurtis", False, True, "annot_"),  # Bray-Curtis with annotation
    ],
)
def test_plot_beta_diversity(
    script_runner,
    test_data_dir,
    amplicon_type,
    level,
    distance_method,
    log10,
    annotated,
    suffix,
):
    """
    Test script for plotting two kinds of beta diversity (with or without annotation)
    """

    # Define test variables
    output_unoise3_dir = "output_unoise3_16s"
    fig_dir = test_data_dir / "figs"

    # Prepare input and output directories
    input_file = test_data_dir / output_unoise3_dir / "unoise3_zotu.biom"
    metadata_file = test_data_dir / "metadata" / f"bulk_{amplicon_type}.tsv"
    taxonomy_file = (
        test_data_dir / output_unoise3_dir / "unoise3_zotu_rrndb_processed.tsv"
    )
    output_dir = fig_dir / f"beta_diversity_{suffix}{amplicon_type}"

    # Construct the base command
    cmd = [
        "plot_dm.py",
        "-i",
        str(input_file),
        "-m",
        str(metadata_file),
        "-t",
        str(taxonomy_file),
        "-f",
        output_dir,
        "-l",
        "otu",
        level,
        "-d",
        distance_method,
        "--hue",
        "rep_group",
    ]

    # Add log10 and transformation if applicable
    if log10:
        cmd.extend(["--log10", "-tr", "relative"])

    # Add annotation option if applicable
    if annotated:
        cmd.append("-a")

    # Run the script
    ret = script_runner.run(cmd)

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"
