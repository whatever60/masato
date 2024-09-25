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
python3 scripts/analysis.py \
  --count_tsv_path outputs/zotu_otus.tsv \
  --sample_metadata_path metadatas/all.tsv \
  --otu_taxonomy_path outputs_rrndb/rrndb_zotus_relabel_processed.tsv \
  --otus_spikein Otu8 Otu96 Otu235 \
  --save_dir outputs_rrndb/all/thres_1 \
  --group_key Group \
  --sample_count_thres 100 \
  --tax_level genus \
  --rel_ab_thres 0.01 \
  --keep_others \
  --keep_unknown
#   --genus_rel_ab_thres 0.05 \
#   --family_rel_ab_thres 0.05
#   --genus_key genus \
#   --family_key family \
"""

import os
import argparse

import pandas as pd
import polars as pl
from skbio.diversity.alpha import shannon


def main():
    parser = argparse.ArgumentParser()

    # Arguments
    parser.add_argument("--count_tsv_path", required=True, help="Path to OTU count tsv.")
    parser.add_argument(
        "--sample_metadata_path", required=True, help="Path to sample metadata."
    )
    parser.add_argument(
        "--otu_taxonomy_path", required=True, help="Path to OTU taxonomy tsv."
    )
    parser.add_argument("--otus_spikein", nargs="+", help="OTUs spikein.", default=[])
    parser.add_argument("--sample_weight_key", default=None)
    parser.add_argument("--save_dir", required=True, help="Directory to save the outputs.")
    parser.add_argument("--group_key", default=None)
    parser.add_argument(
        "--sample_count_thres", type=int, required=True, help="Sample count threshold."
    )
    parser.add_argument(
        "--tax_level",
        type=str,
        required=True,
        help="Taxonomy level of interest, must be a column in the taxonomy table.",
    )
    parser.add_argument(
        "--rel_ab_thres",
        type=float,
        required=True,
        help="Genus relative abundance threshold.",
    )
    # parser.add_argument("--genus_key", required=True, help="Key for genus.")
    # parser.add_argument("--family_key", required=True, help="Key for family.")
    # parser.add_argument(
    #     "--genus_rel_ab_thres",
    #     type=float,
    #     required=True,
    #     help="Genus relative abundance threshold.",
    # )
    # parser.add_argument(
    #     "--family_rel_ab_thres",
    #     type=float,
    #     required=True,
    #     help="Family relative abundance threshold.",
    # )
    parser.add_argument("--keep_others", action="store_true", default=False)
    parser.add_argument("--keep_unknown", action="store_true", default=False)
    args = parser.parse_args()

    # Extract arg values
    count_tsv_path = args.count_tsv_path
    sample_metadata_path = args.sample_metadata_path
    otu_taxonomy_path = args.otu_taxonomy_path
    otus_spikein = args.otus_spikein
    sample_weight_key = args.sample_weight_key
    save_dir = args.save_dir
    # genus_key = args.genus_key
    # family_key = args.family_key
    group_key = args.group_key
    sample_count_thres = args.sample_count_thres
    tax_level = args.tax_level
    rel_ab_thres = args.rel_ab_thres
    # genus_rel_ab_thres = args.genus_rel_ab_thres
    # family_rel_ab_thres = args.family_rel_ab_thres
    keep_others = args.keep_others
    keep_unknown = args.keep_unknown

    os.makedirs(save_dir, exist_ok=True)

    # make sure the first column is sample id
    sample_metadata = pd.read_table(sample_metadata_path, index_col=0)
    if sample_weight_key is None:
        sample_weight_key = "sample_weight"
        sample_metadata[sample_weight_key] = 1
    if group_key is None:
        sample_metadata[group_key] = sample_metadata.index


    otu_taxonomy = pd.read_table(otu_taxonomy_path, index_col="otu")
    # reorder data_count columns by sample metadata
    # num_samples x num_otus
    # exclude OTUs without taxonomy assignment.
    data_count = pd.read_table(count_tsv_path, index_col=0)[sample_metadata.index].T
    common_otus = list(set(data_count.columns).intersection(otu_taxonomy.index))
    data_count = data_count[common_otus]
    otu_taxonomy = otu_taxonomy.loc[common_otus]
    pass_qc = data_count.sum(axis=1) >= sample_count_thres
    sample_metadata["total_reads"] = data_count.drop(otus_spikein, axis=1).sum(axis=1)
    sample_metadata["spikein_reads"] = data_count.loc[:, otus_spikein].sum(axis=1)

    # filter poor samples
    sample_metadata = sample_metadata.loc[pass_qc]
    data_count = data_count.loc[pass_qc].drop(otus_spikein, axis=1)

    if not otus_spikein:
        data_abs = data_count.div(sample_metadata[sample_weight_key], axis=0)
    else:
        data_abs = data_count.div(
            sample_metadata.spikein_reads * sample_metadata[sample_weight_key], axis=0
        )

    # for relative abundance, spike-in and sample weight do not matter.
    data_rel = data_count.div(sample_metadata.total_reads, axis=0)
    sample_metadata["shannon_otu"] = data_rel.apply(shannon, axis=1)
    sample_metadata["richness_otu"] = (data_rel > 0).sum(axis=1)
    data_abs.to_csv(f"{save_dir}/otu_ab_abs.csv")
    data_rel.to_csv(f"{save_dir}/otu_ab_rel.csv")

    melt = data_rel.melt(ignore_index=False).reset_index()
    melt.columns = ["sample", "otu", "rel_ab"]
    melt = melt.merge(sample_metadata, left_on="sample", right_index=True).merge(
        otu_taxonomy, left_on="otu", right_index=True
    )

    melt_mean = (
        pl.from_pandas(melt)
        .groupby(["otu", tax_level, group_key])
        .agg(rel_ab=pl.col("rel_ab").mean())
    )

    tax_keep = (
        melt_mean.groupby(tax_level)
        .agg(ab=pl.max("rel_ab"))
        .filter(pl.col("ab") > rel_ab_thres)
    )[tax_level]

    tax_table_dict = (
        melt_mean.groupby([tax_level, group_key])
        .agg(rel_ab=pl.col("rel_ab").sum())
        .with_columns(is_good=pl.col(tax_level).is_in(tax_keep))
        .partition_by("is_good", as_dict=True)
    )

    tax_table_good = tax_table_dict[True].drop("is_good")
    try:
        tax_table_bad = (
            tax_table_dict[False]
            .with_columns(pl.lit("others").alias(tax_level))
            .drop("is_good")
        )
    except KeyError:
        tax_table = tax_table_good
    else:
        if keep_others:
            tax_table = pl.concat([tax_table_good, tax_table_bad])
        else:
            tax_table = tax_table_good

    if not keep_unknown:
        tax_table = tax_table.filter(pl.col(tax_level) != "unknown")

    tax_table = (
        tax_table.pivot(
            index=group_key,
            columns=tax_level,
            values="rel_ab",
            aggregate_function="sum",
        )
        .to_pandas()
        .set_index(group_key)
    )

    if not keep_others or not keep_unknown:
        # each row no long adds up to one since low abundance OTUs are filtered
        # genus_table = genus_table.div(genus_table.sum(axis=1), axis=0)
        # family_table = family_table.div(family_table.sum(axis=1), axis=0)
        tax_table = tax_table.div(tax_table.sum(axis=1), axis=0)

    # genus_table.to_csv(f"{save_dir}/genus_rel_ab.csv")
    # family_table.to_csv(f"{save_dir}/family_rel_ab.csv")
    tax_table = tax_table.loc[sample_metadata[group_key]]
    tax_table.to_csv(f"{save_dir}/{tax_level}_rel_ab.csv")

    sample_metadata.to_csv(f"{save_dir}/sample_metadata.csv")
