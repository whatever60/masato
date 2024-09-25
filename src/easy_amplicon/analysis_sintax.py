"""
Example run:
python3 analysis_sintax.py \
    --input /home/ubuntu/data/carlotta_troubleshotting/sintax.tsv \
    --metadata /home/ubuntu/data/carlotta_troubleshotting/sample_metadata.tsv \
    --output_dir /home/ubuntu/data/carlotta_troubleshotting/res_sintax \
    --genus_spikein Sporosarcina \
    --group_key Group \
    --sample_count_thres 500 \
    --genus_rel_ab_thres 0.05 \
    --family_rel_ab_thres 0.05
"""

import argparse
from typing import Sequence

import numpy as np
import pandas as pd
import polars as pl
from skbio.diversity.alpha import shannon


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str)
    parser.add_argument("--metadata", type=str)
    parser.add_argument("--output_dir", type=str)
    parser.add_argument("--genus_spikein", nargs="+", required=True, help="OTUs spikein.")
    parser.add_argument("--sample_weight_key", default=None)
    parser.add_argument("--group_key", default=None)
    parser.add_argument(
        "--sample_count_thres", type=int, required=True, help="Sample count threshold."
    )
    parser.add_argument(
        "--genus_rel_ab_thres",
        type=float,
        required=True,
        help="Genus relative abundance threshold.",
    )
    parser.add_argument(
        "--family_rel_ab_thres",
        type=float,
        required=True,
        help="Family relative abundance threshold.",
    )

    args = parser.parse_args()
    input_path = args.input
    metadata_path = args.metadata
    output_dir = args.output_dir
    genus_spikein = args.genus_spikein
    sample_weight_key = args.sample_weight_key
    group_key = args.group_key
    sample_count_thres = args.sample_count_thres
    genus_rel_ab_thres = args.genus_rel_ab_thres
    family_rel_ab_thres = args.family_rel_ab_thres

    sample_metadata = pd.read_table(metadata_path).rename(columns={"Sample": "sample"})
    if sample_weight_key is None:
        sample_weight_key = "sample_weight"
        sample_metadata[sample_weight_key] = 1

    sample_metadata = pl.from_pandas(sample_metadata).sort("sample")

    df = pl.read_csv(
        input_path,
        separator="\t",
        new_columns=["read_id", "prob", "strand", "assign"],
        has_header=False,
    )

    # numpy array is not a Sequence
    Sequence.register(np.ndarray)
    columns = np.char.array(
        ["domain", "phylum", "class", "order", "family", "genus", "species"]
    )

    df_level = df.select(
        pl.col("prob")
        .str.extract_all(r":([^()]+)")
        .list.eval(pl.element().str.lstrip(":"))
        .list.to_struct(fields=columns)
        .alias("level"),
    ).unnest("level")
    df_p = df.select(
        pl.col("prob")
        .str.extract_all(r"\(([\d.]+)\)")
        .list.eval(pl.element().str.slice(1, 4).cast(pl.Float64))
        .list.to_struct(fields=columns + "_p")
        .alias("p"),
    ).unnest("p")
    df_assign = df.select(
        pl.col("assign")
        .str.extract_all(r":([^,]+)")
        .list.eval(pl.element().str.lstrip(":"))
        .list.to_struct(fields=columns + "_conf"),
    ).unnest("assign")

    df_expand = pl.concat(
        [
            df.select(
                "read_id", pl.col("read_id").str.split(".").list.first().alias("sample")
            ),
            df_level,
            df_p,
            df_assign,
        ],
        how="horizontal",
    ).join(sample_metadata, on="sample")

    sample_count_melt = df_expand.groupby(["sample", "genus", "family"]).count()
    sample_qc = (
        sample_count_melt.groupby("sample")
        .agg(pl.sum("count"))
        .filter(pl.col("count") >= sample_count_thres)
        .select("sample")
        .to_series()
    )
    sample_count_melt = sample_count_melt.filter(pl.col("sample").is_in(sample_qc))
    sample_metadata = sample_metadata.filter(pl.col("sample").is_in(sample_qc))

    sample_count_spike_in = (
        sample_count_melt.filter(pl.col("genus").is_in(genus_spikein))
        .groupby("sample")
        .agg(pl.sum("count"))
        .sort("sample")
    )
    sample_count = sample_count_melt.filter(~pl.col("genus").is_in(genus_spikein))
    normalize_factor = sample_count_spike_in["count"] * sample_metadata[sample_weight_key]

    sample_genus_count = (
        sample_count.groupby(["sample", "genus"])
        .agg(pl.sum("count"))
        .pivot(values="count", index="sample", columns="genus")
        .sort("sample")
        .fill_null(0)
        .to_pandas()
        .set_index("sample")
    )
    sample_family_count = (
        sample_count.groupby(["sample", "family"])
        .agg(pl.sum("count"))
        .pivot(values="count", index="sample", columns="family")
        .sort("sample")
        .fill_null(0)
        .to_pandas()
        .set_index("sample")
    )

    sample_metadata = sample_metadata.to_pandas().set_index("sample")
    sample_metadata["spikein_reads"] = sample_count_spike_in["count"]
    sample_metadata["normalize_factor"] = normalize_factor

    # genus/family-level qc
    sample_genus_ab_rel = sample_genus_count.div(sample_genus_count.sum(axis=1), axis=0)
    sample_family_ab_rel = sample_family_count.div(sample_family_count.sum(axis=1), axis=0)
    genus_keep = sample_genus_ab_rel.sum(axis=0) >= genus_rel_ab_thres
    family_keep = sample_family_ab_rel.sum(axis=0) >= family_rel_ab_thres
    sample_genus_count = sample_genus_count.loc[:, genus_keep]
    sample_family_count = sample_family_count.loc[:, family_keep]

    # recalculate genus/family-level abundance
    sample_genus_ab_rel = sample_genus_count.div(sample_genus_count.sum(axis=1), axis=0)
    sample_family_ab_rel = sample_family_count.div(sample_family_count.sum(axis=1), axis=0)
    sample_genus_ab_abs = sample_genus_count.div(normalize_factor, axis=0)
    sample_family_ab_abs = sample_family_count.div(normalize_factor, axis=0)
    group_genus_ab_abs = (
        sample_genus_ab_abs.join(sample_metadata[group_key]).groupby(group_key).mean()
    )
    group_family_ab_abs = (
        sample_family_ab_abs.join(sample_metadata[group_key]).groupby(group_key).mean()
    )
    group_genus_ab_rel = (
        sample_genus_ab_rel.join(sample_metadata[group_key]).groupby(group_key).mean()
    )
    group_family_ab_rel = (
        sample_family_ab_rel.join(sample_metadata[group_key]).groupby(group_key).mean()
    )

    sample_genus_ab_abs.to_csv(f"{output_dir}/sample_genus_abs_ab.csv")
    sample_family_ab_abs.to_csv(f"{output_dir}/sample_family_abs_ab.csv")
    sample_genus_ab_rel.to_csv(f"{output_dir}/sample_genus_rel_ab.csv")
    sample_family_ab_rel.to_csv(f"{output_dir}/sample_family_rel_ab.csv")
    group_genus_ab_abs.to_csv(f"{output_dir}/group_genus_abs_ab.csv")
    group_family_ab_abs.to_csv(f"{output_dir}/group_family_abs_ab.csv")
    group_genus_ab_rel.to_csv(f"{output_dir}/group_genus_rel_ab.csv")
    group_family_ab_rel.to_csv(f"{output_dir}/group_family_rel_ab.csv")

    sample_metadata["shannon_genus"] = sample_genus_ab_rel.apply(shannon, axis=0)
    sample_metadata["shannon_family"] = sample_family_ab_rel.apply(shannon, axis=0)
    sample_metadata.to_csv(f"{output_dir}/sample_metadata.csv")
