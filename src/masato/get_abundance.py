#!/usr/bin/env python3
"""
Take a OTU count table, sample metadata, and OTU taxonomy table as input, output 
tables of relative/absolute abundance at OTU level and relative abundance at 
genus/family level.

QC is also performed during the process. Poor samples with low total OTU counts are 
removed, and poor genus and family with low total relative abundance are removed.

Spike-in genus needs to be specified as an argument. Sample weight could also be 
provided as an argument. Spike-in genus won't be considered in the final tables. 

Absolute OTU abundance is calculated by dividing OTU counts by spike-in counts and 
sample weight.

Relative OTU abundance is calculated by normalizing absolute OTU abundance to 1.

Relative genus/family abundance is calculated by first averaging over any replicates 
within group and averaging again relative OTU abundance of OTUs corresponding to the 
same genus/family.

Note that rows of OTU abundance tables are sample, while rows of genus/family abundance 
tables are groups.

Example run:
get_abundance.py \
  --otu_count_tsv_path outputs_unoise3/unoise3_zotu.tsv \
  --metadata_path metadatas/all.tsv \
  --otu_taxonomy_path outputs_unoise3/unoise3_zotu_rrndb_processed.tsv \
  --output_dir outputs_rrndb/all \
  --spikein_taxa_key spike_in \
  --rep_group_key rep_group \
  --tax_levels otu genus family order \
  --rel_ab_thres 0.01 \
  --keep_others \
  --keep_unknown
"""

import os
import argparse

import pandas as pd
import polars as pl

from masato.utils import read_table, read_table


def _find_otus_by_taxon(df_tax: pd.DataFrame, taxon: str) -> list[str]:
    # index of df_tax are OTU numbers
    level, name = taxon.split(";")
    return df_tax.query(f"{level} == '{name}'").index.tolist()


def _calc_norm_factor(
    df_otu_count: pd.DataFrame,
    df_meta: pd.DataFrame,
    df_tax: pd.DataFrame,
    sample_weight_key: str = None,
    spikein_taxa_key: str = None,
) -> pd.DataFrame:
    """
    Get absolute and relative abundance of OTUs. In this step replicate samples are
        still not averaged.

    The index of `df_otu_count` should be sample names, and columns should be OTU
        numbers.
    The index of `df_tax` should be OTU numbers.
    The index of `df_meta` should be sample names. Only samples specified in `df_meta`
        will be kept in the output.

    A new metadata dataframe will be returned.
    """
    spikein2otus = find_spikein_otus(df_meta, df_tax, spikein_taxa_key)

    def compute_spikein_and_norm(row: pd.Series) -> pd.Series:
        otu_count_series = df_otu_count.loc[row.name]

        if spikein_taxa_key and isinstance(row[spikein_taxa_key], str):
            otus_spikein = [
                otu
                for spikein in row[spikein_taxa_key].split(",")
                for otu in spikein2otus[spikein]
            ]
        else:
            otus_spikein = []
        if otus_spikein:
            spikein_reads = otu_count_series.loc[otus_spikein].sum()
            otu_count_series = otu_count_series.drop(otus_spikein)
        else:
            spikein_reads = -1

        norm_factor = spikein_reads if otus_spikein else 1
        if sample_weight_key and row[sample_weight_key]:
            norm_factor *= row[sample_weight_key]
        return pd.Series(
            {
                "spikein_reads": spikein_reads,
                "non_spikein_reads": otu_count_series.sum(),
                "norm_factor": norm_factor,
            }
        )

    if df_meta.shape[1]:
        df_meta_add = df_meta.apply(compute_spikein_and_norm, axis=1)
    else:
        df_meta_add = pd.DataFrame(
            {
                "spikein_reads": -1,
                "non_spikein_reads": df_otu_count.sum(axis=1),
                "norm_factor": 1,
            }
        )
    df_meta = df_meta.join(df_meta_add)
    return df_meta


# def _calc_alpha_metrics(df: pd.DataFrame) -> pd.DataFrame:
#     """Calculate alpha diversity metrics for each sample.
#     Note that you should ignore chao1 if your input data is not integer count.
#     """
#     return pd.DataFrame(
#         {
#             "chao1": df.apply(chao1, axis=1),
#             "richness": df.apply(lambda x: (x > 0).sum(), axis=1),
#             "shannon": df.apply(shannon, axis=1),
#             "simpson": df.apply(simpson, axis=1),
#         },
#         index=df.index,
#     )


def _agg_along_axis(
    df: pd.DataFrame, series: pd.Series, axis: int, aggfunc: str = None
) -> pd.DataFrame:
    """
    Aggregate a DataFrame along a specified axis using a given function.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame to be aggregated.
    series : pd.Series
        A Series used for grouping the DataFrame.
    axis : int
        The axis along which to perform the aggregation.
        0 for rows and 1 for columns.
    aggfunc : str, optional
        The aggregation function to use. Must be either 'mean' or 'sum'.
        If None, the function defaults to mean for axis 0 and sum for axis 1.

    Returns
    -------
    pd.DataFrame
        The aggregated DataFrame with groups based on the series provided.

    Raises
    ------
    ValueError
        If `aggfunc` is not 'mean' or 'sum'.
    """
    index_name = df.index.name or "index"
    series_index_name = series.index.name or "index"
    df = pl.from_pandas(df.reset_index().rename(columns={index_name: "row"}))
    exprs = [pl.col("rel_ab").mean(), pl.col("rel_ab").sum()]

    if not aggfunc:
        expr = exprs[axis]
    elif aggfunc == "sum":
        expr = exprs[1]
    elif aggfunc == "mean":
        expr = exprs[0]
    else:
        raise ValueError("expr must be either 'mean' or 'sum'.")

    if axis == 0:
        series = pl.from_pandas(
            series.reset_index(name="group").rename(columns={series_index_name: "row"})
        )
        df_g = (
            df.unpivot(index="row", variable_name="col", value_name="rel_ab")
            .join(series, on="row")
            .group_by(["col", "group"])
            .agg(rel_ab=expr)
            .pivot(
                index="group",
                on="col",
                values="rel_ab",
            )
            .to_pandas()
            .set_index("group")
        )
    else:
        series = pl.from_pandas(
            series.reset_index(name="group").rename(columns={series_index_name: "col"})
        )
        df_g = (
            df.unpivot(index="row", variable_name="col", value_name="rel_ab")
            .join(series, on="col")
            .group_by(["row", "group"])
            .agg(rel_ab=expr)
            .pivot(
                index="group",
                on="row",
                values="rel_ab",
            )
            .to_pandas()
            .set_index("group")
        ).T
    # on, groupby = ("row", "col") if axis == 0 else ("col", "row")
    # df_g = (
    #     df.melt(id_vars="row", variable_name="col", value_name="rel_ab")
    #     .join(series, on=on)
    #     .group_by([groupby, "group"])
    #     .agg(rel_ab=pl.col("rel_ab").mean())
    #     .pivot(
    #         index="group",
    #         columns=groupby,
    #         values="rel_ab",
    #     )
    #     .to_pandas()
    #     .set_index("group")
    # )
    # df_g = df_g if axis == 0 else df_g.T
    df_g.index.name = index_name
    return df_g


def _taxa_qc(
    df_tax_rel_ab: pd.DataFrame,
    rel_ab_thres: float,
    keep_rare: bool,
    keep_unknown: bool,
) -> pd.DataFrame:
    """
    Perform quality control on taxonomic data by filtering rare and unknown taxa.

    Parameters
    ----------
    df_tax_rel_ab : pd.DataFrame
        A DataFrame containing the relative abundance of taxa. Columns represent
            different taxa, and rows represent samples.
    rel_ab_thres : float
        Threshold to determine rare taxa. If float, taxa with relative abundance below
            this threshold in all samples are considered rare.
        If int, taxa with a rank of relative abundance above this threshold in all
            samples are considered rare.
    keep_rare : bool
        If False, remove rare taxa. If True, combine rare taxa into a single 'others' column.
    keep_unknown : bool
        If False, remove taxa labeled as unknown. If True, keep unknown taxa.

    Returns
    -------
    pd.DataFrame
        A DataFrame with filtered taxa based on the provided thresholds and options.

    Notes
    -----
    A taxon is defined as rare:
    - When `rel_ab_thres` is a float: If its relative abundance is below the threshold in all samples.
    - When `rel_ab_thres` is an int: If the rank of its relative abundance is above the threshold in all samples.

    A taxon is defined as unknown if its name contains 'unknown' or ends with '_UNKNOWN'.
    """
    df_tax_rel_ab_temp = df_tax_rel_ab[
        [
            t
            for t in df_tax_rel_ab.columns
            if not t.endswith("_UNKNOWN") and t != "unknown"
        ]
    ]
    if rel_ab_thres and int(rel_ab_thres) == rel_ab_thres:
        rare_taxa = df_tax_rel_ab_temp.columns[
            (
                df_tax_rel_ab_temp.rank(axis=1, ascending=False, method="max")
                > rel_ab_thres
            ).all(axis=0)
        ].tolist()
    else:
        rare_taxa = df_tax_rel_ab_temp.columns[
            (df_tax_rel_ab_temp < rel_ab_thres).all(axis=0)
        ].tolist()
    # if "unknown" in rare_taxa:
    #     rare_taxa.remove("unknown")
    # rare_taxa = [t for t in rare_taxa if not t.endswith("_UNKNOWN") or t == "unknown"]
    if not keep_rare:
        df_tax_rel_ab = df_tax_rel_ab.drop(rare_taxa, axis=1)
    elif rare_taxa:
        df_tax_rel_ab_rare = (
            df_tax_rel_ab[rare_taxa].sum(axis=1).to_frame(name="others")
        )
        df_tax_rel_ab = df_tax_rel_ab.drop(rare_taxa, axis=1)
        df_tax_rel_ab = pd.concat([df_tax_rel_ab, df_tax_rel_ab_rare], axis=1)
    if not keep_unknown:
        df_tax_rel_ab = df_tax_rel_ab[
            df_tax_rel_ab.columns[
                ~df_tax_rel_ab.columns.str.lower().str.endswith("unknown")
            ]
        ]
        if (df_tax_rel_ab.dtypes != "int").any():
            df_tax_rel_ab = df_tax_rel_ab.div(df_tax_rel_ab.sum(1), 0).fillna(0)
    return df_tax_rel_ab


def read_tables(
    otu_count_table: str| pd.DataFrame,
    metadata_path: str | pd.DataFrame| None = None,
    otu_taxonomy_path: str| pd.DataFrame| None = None,
    warning: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    This function reads several tables from amplicon analysis, does necessary
        preprocessing, and return tables as dataframes.

    It reads OTU count table output by `VSEARCH` or `USERACH` where rows are sample
        names and columns are OTU names. The header of the first column is `#OTU ID`.

    Optionally, it takes a metadata table whose first column should be sample name. It
        ensures each sample in the metadata has associated OTU counts. If some samples
        are absent in the OTU count table, it fills in missing samples with all zero
        counts. Samples not in the metadata will be removed from the OTU count table.
        Samples in the returned OTU count table will be reordered according to the
        metadata.

    Optionally, it takes a taxonomy table whose the first column is OTU names and
        columns are taxonomy levels, and each taxonomy level has a paired confidence
        column. For example, the level `geuns` has a paired confidence column `geuns_p`.
    This function takes the intersection of OTUs in the OTU count table and the taxonomy
        table. If some OTUs do not have taxonomy, their count will be aggregated to a
        new "unknown" column, and a new row will be added to the taxonomy dataframe
        where taxonomy levels are "unknown" and confidences are -1.
    This function also copies index to a new column `otu` and `otu_p` if not already
        present.

    If sample metadata is not provided, metadata dataframe will be initialized from the
        samples in the OTU count table and have no columns.
    If taxonomy is not provided, taxonomy dataframe will be initialized from the OTUs in
        the OTU count table, and with only `otu` and `otu_p` columns.

    Returns:
        OTU count table: pd.DataFrame. Rows are samples, columns are OTUs, transposed
            from the input table.
        Sample metadata: pd.DataFrame.
        OTU taxonomy: pd.DataFrame.
    """
    if isinstance(otu_count_table, str):
        df_otu_count = read_table(otu_count_table, index_col="#OTU ID")
    else:
        df_otu_count = otu_count_table.copy()
    if not df_otu_count.columns.is_unique:
        raise ValueError("Sample names in OTU count table must be unique.")

    if metadata_path is not None:
        if isinstance(metadata_path, str):
            df_meta = read_table(metadata_path, index_col="sample", comment="#")
        else:
            df_meta = metadata_path.copy()
        if not df_meta.index.is_unique:
            raise ValueError("Sample names in metadata must be unique.")

        # add missing samples with all zero counts
        missing_samples = list(set(df_meta.index) - set(df_otu_count.columns))
        if missing_samples:
            if warning:
                print(
                    f"WARNING: The following {len(missing_samples)} samples are "
                    "missing from OTU count table:"
                )
                print("\t", end="")
                print(", ".join(map(str, sorted(missing_samples))))
                print("Filling in missing samples with all zero counts.")
            df_otu_count = pd.concat(
                [
                    df_otu_count,
                    pd.DataFrame(0, index=df_otu_count.index, columns=missing_samples),
                ],
                axis=1,
            )
        # only care about samples in metadata
        df_otu_count = df_otu_count[df_meta.index].transpose()
    else:
        df_otu_count = df_otu_count.transpose()
        df_meta = pd.DataFrame(index=df_otu_count.index)

    if otu_taxonomy_path is not None:
        # df_tax = pd.read_table(otu_taxonomy_path, index_col="otu").rename(
        #     {"otu.1": "otu"}, axis=1
        # )
        if isinstance(otu_taxonomy_path, str):
            df_tax = read_table(otu_taxonomy_path, index_col="otu")
        else:
            df_tax = otu_taxonomy_path.copy()
        if not df_tax.index.is_unique:
            raise ValueError("OTU numbers in taxonomy table must be unique.")
        # only care about OTUs with taxonomy
        otu_in_tax = set(df_tax.index)
        otu_in_count = set(df_otu_count.columns)
        no_tax_otus = list(otu_in_count - otu_in_tax)
        no_count_otus = list(otu_in_tax - otu_in_count)
        if no_tax_otus:
            common_otus = list(otu_in_tax & otu_in_count)
            df_tax = df_tax.loc[common_otus]
            df_tax_add = pd.DataFrame(
                "unknown",
                index=no_tax_otus,
                columns=df_tax.columns,
            )
            df_tax_add.loc[:, df_tax.columns[df_tax.columns.str.endswith("_p")]] = -1
            df_tax = pd.concat([df_tax, df_tax_add]).loc[df_otu_count.columns]
        if no_count_otus:  # add missing OTUs with all zero counts
            df_otu_count = pd.concat(
                [
                    df_otu_count,
                    pd.DataFrame(0, columns=no_count_otus, index=df_otu_count.index),
                ],
                axis=1,
            )
    else:
        df_tax = pd.DataFrame(index=df_otu_count.columns)
    if "otu" not in df_tax.columns:
        df_tax["otu"] = df_tax.index
    if "otu_p" not in df_tax.columns:
        df_tax["otu_p"] = 1
    return df_otu_count, df_meta, df_tax


def find_spikein_otus(
    df_meta: pd.DataFrame,
    df_tax: pd.DataFrame,
    spikein_taxa_key: str = None,
    only_otus: bool = False,
) -> dict[str, list[str]] | list[str]:
    if spikein_taxa_key is None:
        return []
    try:
        all_spikeins = (
            df_meta[spikein_taxa_key]
            .dropna()
            .str.split(",")
            .explode()
            .unique()
            .tolist()
        )
    except AttributeError:  # happens when all spikein_taxa_key are NaN
        all_spikeins = []
    spikein2otus = {
        spikein: _find_otus_by_taxon(df_tax, spikein) for spikein in all_spikeins
    }
    if only_otus:
        return list(set([otu for otus in spikein2otus.values() for otu in otus]))
    else:
        return spikein2otus


def get_otu_count(
    otu_count_table: str | pd.DataFrame,
    metadata_path: str | pd.DataFrame | None = None,
    otu_taxonomy_path: str | pd.DataFrame | None = None,
    sample_weight_key: str | None = None,
    spikein_taxa_key: str | None = None,
    warning: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df_otu_count, df_meta, df_tax = read_tables(
        otu_count_table, metadata_path, otu_taxonomy_path, warning=warning
    )
    assert len(df_otu_count) == len(df_meta) and len(df_otu_count.columns) == len(df_tax)
    df_meta = _calc_norm_factor(
        df_otu_count,
        df_meta,
        df_tax,
        sample_weight_key,
        spikein_taxa_key,
    )
    df_meta["sequencing_depth"] = df_otu_count.sum(
        axis=1
    ) - df_meta.spikein_reads.to_numpy().clip(0)
    spikein_otus = find_spikein_otus(df_meta, df_tax, spikein_taxa_key, only_otus=True)
    df_otu_count = df_otu_count.drop(spikein_otus, axis=1)
    df_tax = df_tax.drop(spikein_otus, axis=0)  
    return df_otu_count, df_meta, df_tax


def main():
# if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Arguments
    parser.add_argument(
        "-i", "--otu_count_tsv_path", required=True, help="Path to OTU count tsv."
    )
    parser.add_argument(
        "-m", "--metadata_path", required=True, help="Path to sample metadata."
    )
    parser.add_argument(
        "-t", "--otu_taxonomy_path", required=True, help="Path to OTU taxonomy tsv."
    )
    parser.add_argument("-sp", "--spikein_taxa_key", type=str, default="spike_in")
    parser.add_argument("-r", "--rep_group_key", type=str, default="rep_group")
    parser.add_argument("-w", "--sample_weight_key", default=None, type=str)
    parser.add_argument(
        "-o", "--output_dir", required=True, help="Directory to save the outputs."
    )
    parser.add_argument(
        "-l",
        "--tax_levels",
        nargs="+",
        required=True,
        help="Taxonomy level of interest, must be a column in the taxonomy table.",
    )
    parser.add_argument(
        "-a",
        "--rel_ab_thresholds",
        nargs="+",
        type=float,
        help="Genus relative abundance threshold.",
    )
    parser.add_argument("--keep_others", action="store_true", default=False)
    parser.add_argument("--keep_unknown", action="store_true", default=False)
    args = parser.parse_args()

    # Extract arg values
    count_tsv_path = args.otu_count_tsv_path
    metadata_path = args.metadata_path
    otu_taxonomy_path = args.otu_taxonomy_path
    spikein_taxa_key = args.spikein_taxa_key
    sample_weight_key = args.sample_weight_key
    output_dir = args.output_dir
    rep_group_key = args.rep_group_key
    tax_levels = args.tax_levels
    rel_ab_thresholds = args.rel_ab_thresholds
    keep_others = args.keep_others
    keep_unknown = args.keep_unknown

    os.makedirs(output_dir, exist_ok=True)

    # read inputs
    df_otu_count, df_meta, df_tax = get_otu_count(
        count_tsv_path,
        metadata_path,
        otu_taxonomy_path,
        sample_weight_key,
        spikein_taxa_key,
    )

    df_meta.to_csv(f"{output_dir}/metadata_sample.csv")
    df_otu_abs_ab = df_otu_count.div(df_meta.norm_factor, axis=0)
    df_otu_abs_ab.to_csv(f"{output_dir}/ab_abs_s_otu.csv")

    df_otu_rel_ab = df_otu_count.div(df_otu_count.sum(axis=1), axis=0)
    df_otu_rel_ab_g = _agg_along_axis(df_otu_rel_ab, df_meta[rep_group_key], axis=0)

    # df_otu_count.to_csv(f"{output_dir}/count_sample_otu.csv")
    # absolute abundance is only for OTU level

    if len(rel_ab_thresholds) == 1:
        rel_ab_thresholds = rel_ab_thresholds * len(tax_levels)
    for level, rel_ab_thres in zip(tax_levels, rel_ab_thresholds):
        df_tax_rel_ab_g = _taxa_qc(
            _agg_along_axis(df_otu_rel_ab_g, df_tax[level], axis=1),
            rel_ab_thres,
            keep_others,
            keep_unknown,
        )
        df_tax_rel_ab = _taxa_qc(
            _agg_along_axis(df_otu_rel_ab, df_tax[level], axis=1),
            rel_ab_thres,
            keep_others,
            keep_unknown,
        )
        df_tax_rel_ab_g.to_csv(f"{output_dir}/ab_rel_g_{level}.csv")
        df_tax_rel_ab.to_csv(f"{output_dir}/ab_rel_s_{level}.csv")
