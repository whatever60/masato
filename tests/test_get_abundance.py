import pandas as pd
import pytest

from masato.get_abundance import find_spikein_otus, get_otu_count
from masato.isolate_utils import read_isolate_metadata_rich


@pytest.mark.parametrize(
    "taxon, expected",
    [
        ("class;Chloroplast", ["plastid"]),
        ("genus;O'Brien", ["quoted"]),
        ("genus;Bacillus", ["bacterium"]),
    ],
)
def test_find_spikein_otus_handles_reserved_ranks_and_literal_names(taxon, expected):
    taxonomy = pd.DataFrame(
        {
            "class": ["Chloroplast", "Bacilli", "Bacilli"],
            "genus": ["unclassified", "O'Brien", "Bacillus"],
        },
        index=["plastid", "quoted", "bacterium"],
    )
    metadata = pd.DataFrame({"spike_in": [taxon]}, index=["sample1"])

    result = find_spikein_otus(metadata, taxonomy, "spike_in")

    assert result[taxon] == expected


def test_isolate_metadata_excludes_chloroplast_reads_from_bacterial_counts(tmp_path):
    picking_file = tmp_path / "Destination D001 - picked isolates.csv"
    picking_file.write_text(
        'export metadata\nplate metadata\npicking_coord,src_plate,dest_well\n'
        '"(1.0, 2.0)",SRC_1,A1\n'
    )
    metadata = read_isolate_metadata_rich(str(tmp_path))
    counts = pd.DataFrame(
        {"D001_A1": [100, 80, 20, 10]},
        index=["spike", "plastid", "bacterium1", "bacterium2"],
    )
    taxonomy = pd.DataFrame(
        {
            "genus": ["Sporosarcina", "unclassified", "Bacillus", "O'Brien"],
            "class": ["Bacilli", "Chloroplast", "Bacilli", "Bacilli"],
            "family": ["Planococcaceae", "unclassified", "Bacillaceae", "other"],
        },
        index=counts.index,
    )

    filtered, sample_metadata, filtered_taxonomy = get_otu_count(
        counts, metadata, taxonomy, spikein_taxa_key="spike_in_16s"
    )

    assert filtered.columns.tolist() == ["bacterium1", "bacterium2"]
    assert filtered.loc["D001_A1"].tolist() == [20, 10]
    assert filtered_taxonomy.index.tolist() == ["bacterium1", "bacterium2"]
    assert sample_metadata.loc["D001_A1", "sequencing_depth"] == 30
    assert sample_metadata.loc["D001_A1", "non_spikein_reads"] == 30
