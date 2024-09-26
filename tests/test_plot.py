import pytest
from pathlib import Path


@pytest.fixture
def test_data_dir():
    return Path(__file__).parent / "data"


@pytest.mark.parametrize("amplicon_type", ["16s"])
@pytest.mark.parametrize("level", ["genus"])
@pytest.mark.parametrize(
    "fig_dir_order_pair",
    [("abc", "alphabetical"), ("clust", "hierarchical"), ("tax", "taxonomy_tree")],
)
@pytest.mark.parametrize("threshold_suffix_pair", [(1e-2, "1e-2"), (3, "top3")])
@pytest.mark.parametrize("subcommand", ["heatmap", "heatmap_log10", "stacked_bar"])
def test_plot(
    script_runner,
    test_data_dir,
    fig_dir_order_pair,
    subcommand,
    threshold_suffix_pair,
    amplicon_type,
    level,
):
    """
    Test script that mimics a complex nested loop from the original bash script.
    """

    # Unpack the parameters
    (fig_dir_name, order_option) = fig_dir_order_pair
    (abundance_threshold, suffix) = threshold_suffix_pair

    # Define test variables
    output_unoise3_dir = "output_unoise3_16s"

    # Prepare input and output directories
    input_file = test_data_dir / output_unoise3_dir / "unoise3_zotu.biom"
    metadata_file = test_data_dir / "metadata" / f"bulk_{amplicon_type}.tsv"
    taxonomy_file = (
        test_data_dir / output_unoise3_dir / "unoise3_zotu_rrndb_processed.tsv"
    )
    fig_dir = test_data_dir / "figs"
    adjusted_fig_dir = f"{fig_dir_name}_{suffix}"

    # Construct the command to run
    cmd = [
        "plot.py",
        "abundance_group",
        subcommand,
        "-i",
        str(input_file),
        "-m",
        str(metadata_file),
        "-t",
        str(taxonomy_file),
        "-f",
        str(fig_dir / f"abundance_{adjusted_fig_dir}_{amplicon_type}"),
        "-l",
        level,
        "-fo",
        order_option,
        "-a",
        str(abundance_threshold),
    ]

    # Run the script
    ret = script_runner.run(cmd)

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"


@pytest.mark.parametrize("amplicon_type", ["16s"])
@pytest.mark.parametrize("level", ["genus"])
@pytest.mark.parametrize("rarefying", [False, True])  # To handle both commands
def test_plot_alpha_diversity(script_runner, test_data_dir, amplicon_type, level, rarefying):
    """
    Test script for plotting alpha diversity, with and without rarefying.
    """

    # Define test variables
    output_unoise3_dir = "output_unoise3_16s"
    fig_dir = test_data_dir / "figs"

    # Prepare input and output directories
    input_file = test_data_dir / output_unoise3_dir / "unoise3_zotu.biom"
    metadata_file = test_data_dir / "metadata" / f"bulk_{amplicon_type}.tsv"
    taxonomy_file = test_data_dir / output_unoise3_dir / "unoise3_zotu_rrndb_processed.tsv"

    # Construct the base command
    cmd = [
        "plot.py",
        "stats_sample_count",
        "-i", str(input_file),
        "-m", str(metadata_file),
        "-t", str(taxonomy_file),
        "-f", str(fig_dir / f"alpha_diversity_{'rarefying_' if rarefying else ''}{amplicon_type}"),
        "-l", "otu", level
    ]

    # Add rarefying options if needed
    if rarefying:
        cmd.extend(["-rv", "500", "-rn", "50"])

    # Run the script
    ret = script_runner.run(cmd)

    # Assertions
    assert ret.success, f"Script failed with return code {ret.returncode}"
    # output_file = fig_dir / f"alpha_diversity_{'rarefying_' if rarefying else ''}{amplicon_type}.png"
    # assert output_file.exists(), f"Output file not created: {output_file}"
