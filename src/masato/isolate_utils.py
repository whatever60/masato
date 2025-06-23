#!/usr/bin/env python
import argparse
from collections import defaultdict
import os
import glob

import numpy as np
import pandas as pd

from masato.get_abundance import get_otu_count, _agg_along_axis, _taxa_qc
from masato.utils import write_table


def _read_isolate_metadata(isolate_metadata_dir: str) -> pd.DataFrame:
    isolate_metadata_list = []
    for path in glob.glob(
        os.path.join(isolate_metadata_dir, "Destination * - * *.csv")
    ):
        dest_plate_barcode = os.path.basename(path).split(" ")[1]
        df = pd.read_csv(path, skiprows=2).drop_duplicates()
        df.iloc[:, 1] = df.iloc[:, 1].map(lambda x: x.split("_")[0])
        df["dest_plate"] = dest_plate_barcode
        isolate_metadata_list.append(df)
    return pd.concat(isolate_metadata_list, ignore_index=True)


def read_isolate_metadata_rich(
    isolate_metadata_dir: str, plate_metadata_path: str | None = None
) -> pd.DataFrame:
    isolate_metadata = _read_isolate_metadata(isolate_metadata_dir)
    isolate_metadata.columns = ["picking_coord", "src_plate", "dest_well", "dest_plate"]
    if plate_metadata_path is not None:
        plate_metadata = pd.read_csv(plate_metadata_path)
        if not plate_metadata.barcode.is_unique:
            raise ValueError("Plate metadata barcode column is not unique")
        isolate_metadata = isolate_metadata.merge(
            plate_metadata[
                ["barcode", "medium_type", "sample_type", "plate_name", "sample_group"]
            ].rename({"barcode": "src_plate"}, axis=1),
            on="src_plate",
            how="left",
        )
    else:
        isolate_metadata["medium_type"] = "unknown"
        isolate_metadata["sample_type"] = "unknown"
        isolate_metadata["plate_name"] = "unknown"
        isolate_metadata["sample_group"] = "unknown"
    isolate_metadata["dest_well_barcode"] = (
        isolate_metadata.dest_plate + "_" + isolate_metadata.dest_well
    )
    isolate_metadata["sample"] = isolate_metadata["dest_well_barcode"]
    isolate_metadata[["src_x", "src_y"]] = (
        isolate_metadata["picking_coord"]
        .str.extract(r"\((\d+.\d+), (\d+.\d+)\)")
        .astype(float)
        .to_numpy()
    )
    isolate_metadata["spike_in_16s"] = "genus;Sporosarcina,family;Chloroplast"
    isolate_metadata["spike_in_its"] = None
    return isolate_metadata.sort_values(
        ["sample_group", "sample_type", "medium_type", "src_plate", "sample"]
    ).set_index("sample")



def get_top_cols(
    df: pd.DataFrame,
    top_k: int = 5,
) -> pd.DataFrame:
    ret = {}
    argmax_arr = np.argsort(df, axis=1)
    max_arr = np.sort(df, axis=1)
    for i in range(min(top_k, df.shape[1])):
        ret[f"top_{i + 1}"] = df.columns[argmax_arr[:, -1 - i]]
        ret[f"top_{i + 1}_val"] = max_arr[:, -1 - i]
    return pd.DataFrame(ret, index=df.index)


def _calc_purity(
    df_zotu_count: pd.DataFrame, df_taxonomy: pd.DataFrame, level: str = "otu"
) -> pd.DataFrame:
    """Get purity table for input count table at a given taxonomy level"""
    df_rab = _taxa_qc(
        _agg_along_axis(
            df_zotu_count.div(df_zotu_count.sum(axis=1), axis=0).fillna(0),
            df_taxonomy[level],
            axis=1,
        ),
        rel_ab_thres=0,
        keep_rare=True,
        keep_unknown=True,
    )
    df_purity = get_top_cols(df_rab)
    df_purity["total_counts"] = df_zotu_count.sum(axis=1)
    return df_purity


def get_isolate_purity(
    df_zotu_cb: pd.DataFrame | None = None,
    df_zotu_cf: pd.DataFrame | None = None,
    df_taxonomy_b: pd.DataFrame | None = None,
    df_taxonomy_f: pd.DataFrame | None = None,
    levels: list = ["otu"],
) -> dict[str, list]:
    """Get purity table for multiple input count tables and the combined count table at
    multiple taxonomy levels

    If bacteria zotu count table is given, get its purity.
    If fungi zotu count table is given, get its purity.
    If both are given, get the purity for bacteria and fungi individually and collectively.
    """
    ret_purity = defaultdict(list)
    for level in levels:
        both = 0
        for df_zotu, df_tax in zip(
            [df_zotu_cb, df_zotu_cf], [df_taxonomy_b, df_taxonomy_f]
        ):
            if df_zotu is not None:
                df_purity = _calc_purity(df_zotu, df_tax, level)
                both += 1
                ret_purity[level].append(df_purity)
            else:
                ret_purity[level].append(None)
        if both == 2:
            df_zotu_count = pd.concat([df_zotu_cb, df_zotu_cf], axis=1).fillna(0)
            df_taxonomy = pd.concat([df_taxonomy_b, df_taxonomy_f], axis=0)
            df_purity = _calc_purity(df_zotu_count, df_taxonomy, level)
            ret_purity[level].append(df_purity)
        else:
            ret_purity[level].append(None)
    return ret_purity


# def get_isolate_purity(
#     df_isolate_meta_tsv: str,
#     df_zotu_count_bacteria_tsv: str = None,
#     df_zotu_count_fungi_tsv: str = None,
#     df_taxonomy_bacteria_tsv: str = None,
#     df_taxonomy_fungi_tsv: str = None,
#     level: str = "otu",
# ):
#     """Get purity table for isolates.

#     If bacteria zotu count table is given, get its purity.
#     If fungi zotu count table is given, get its purity.
#     If both are given, get the purity for bacteria and fungi individually and collectively.

#     If taxonomy tables are given, attach the taxonomy information to the purity table.
#     """
#     both = 0
#     ret_count = []
#     ret_tax = []
#     ret_purity = []
#     if df_zotu_count_bacteria_tsv is not None:
#         df_zotu_cb, _, df_taxonomy_b = get_otu_count(
#             df_zotu_count_bacteria_tsv,
#             df_isolate_meta_tsv,
#             df_taxonomy_bacteria_tsv,
#             warning=False,
#         )
#         df_rab_b = _taxa_qc(
#             _agg_along_axis(
#                 df_zotu_cb.div(df_zotu_cb.sum(axis=1), axis=0).fillna(0),
#                 df_taxonomy_b[level],
#                 axis=1,
#             ),
#             rel_ab_thres=0,
#             keep_rare=True,
#             keep_unknown=False,
#         )
#         df_purity_bacteria = get_top_cols(df_rab_b)
#         df_purity_bacteria["total_counts"] = df_zotu_cb.sum(axis=1)
#         both += 1
#         ret_count.append(df_zotu_cb)
#         ret_tax.append(df_taxonomy_b)
#         ret_purity.append(df_purity_bacteria)
#     else:
#         ret_count.append(None)
#         ret_tax.append(None)
#         ret_purity.append(None)
#     if df_zotu_count_fungi_tsv is not None:
#         df_zotu_cf, _, df_taxonomy_f = get_otu_count(
#             df_zotu_count_fungi_tsv,
#             df_isolate_meta_tsv,
#             df_taxonomy_fungi_tsv,
#             warning=False,
#         )
#         df_rab_f = _taxa_qc(
#             _agg_along_axis(
#                 df_zotu_cf.div(df_zotu_cf.sum(axis=1), axis=0).fillna(0),
#                 df_taxonomy_f[level],
#                 axis=1,
#             ),
#             rel_ab_thres=0,
#             keep_rare=True,
#             keep_unknown=False,
#         )
#         df_purity_fungi = get_top_cols(df_rab_f)
#         df_purity_fungi["total_counts"] = df_zotu_cf.sum(axis=1)
#         both += 1
#         ret_count.append(df_zotu_cf)
#         ret_tax.append(df_taxonomy_f)
#         ret_purity.append(df_purity_fungi)
#     else:
#         ret_count.append(None)
#         ret_tax.append(None)
#         ret_purity.append(None)
#     if both == 2:
#         df_zotu_count = pd.concat([df_zotu_cb, df_zotu_cf], axis=1).fillna(0)
#         df_taxonomy = pd.concat([df_taxonomy_b, df_taxonomy_f], axis=0)
#         df_rab = _taxa_qc(
#             _agg_along_axis(
#                 df_zotu_count.div(df_zotu_count.sum(axis=1), axis=0).fillna(0),
#                 df_taxonomy[level],
#                 axis=1,
#             ),
#             rel_ab_thres=0,
#             keep_rare=True,
#             keep_unknown=False,
#         )
#         df_purity = get_top_cols(df_rab)
#         df_purity["total_counts"] = df_zotu_count.sum(axis=1)
#         ret_count.append(df_zotu_count)
#         ret_tax.append(df_taxonomy)
#         ret_purity.append(df_purity)
#     else:
#         ret_count.append(None)
#         ret_tax.append(None)
#         ret_purity.append(None)
#     [df.sort_index(inplace=True) for df in ret_purity if df is not None]
#     try:
#         ret_count = pd.concat([df for df in ret_count if df is not None], axis=1)
#         ret_tax = pd.concat([df for df in ret_tax if df is not None], axis=0)
#     except ValueError:
#         ret_count = ret_tax = None
#     # drop columns of ret_count and rows of ret_tax whose suffix is ZOTU_UNKNWON
#     ret_count = ret_count.loc[:, ~ret_count.columns.str.endswith("ZOTU_UNKNOWN")].copy()
#     ret_tax = ret_tax.loc[~ret_tax.index.str.endswith("ZOTU_UNKNOWN"), :].copy()
#     return ret_count, ret_tax, pd.concat(ret_purity, axis=1)


def combine_count(
    df_isolate_meta_tsv: str,
    df_zotu_count_bacteria_tsv: str = None,
    df_zotu_count_fungi_tsv: str = None,
    df_taxonomy_bacteria_tsv: str = None,
    df_taxonomy_fungi_tsv: str = None,
) -> tuple[list[pd.DataFrame], list[pd.DataFrame]]:
    """This function simply combines the zotu count tables from two different
    experiments. We use it to merge amplicon sequencing data using 16S and ITS primer.

    It preprocesses the data according to what `get_otu_count` does, and
    renames ZOTU names that overlap between the two experiments.
    """
    dfs_zotu_c = []
    dfs_taxonomy = []
    for df_zotu_count_tsv, df_taxonomy_tsv in zip(
        [df_zotu_count_bacteria_tsv, df_zotu_count_fungi_tsv],
        [df_taxonomy_bacteria_tsv, df_taxonomy_fungi_tsv],
    ):
        if df_zotu_count_tsv is not None:
            df_zotu_c, _, df_taxonomy = get_otu_count(
                df_zotu_count_tsv, df_isolate_meta_tsv, df_taxonomy_tsv, warning=False
            )
            dfs_zotu_c.append(df_zotu_c)
            dfs_taxonomy.append(df_taxonomy)
        else:
            dfs_zotu_c.append(None)
            dfs_taxonomy.append(None)

    # rename ZOTU names by prefixing with 16S- for bacteria and ITS- for fungi
    if all(df is not None for df in dfs_zotu_c):
        overlap = set(dfs_zotu_c[0].columns) & set(dfs_zotu_c[1].columns)
        rename_dict_16s = {x: "16S-" + x for x in overlap}
        rename_dict_its = {x: "ITS-" + x for x in overlap}
        dfs_zotu_c[0].rename(columns=rename_dict_16s, inplace=True)
        dfs_zotu_c[1].rename(columns=rename_dict_its, inplace=True)
        dfs_taxonomy[0].rename(index=rename_dict_16s, inplace=True)
        dfs_taxonomy[1].rename(index=rename_dict_its, inplace=True)

    return dfs_zotu_c, dfs_taxonomy


def quality_control(
    df_purity: list[pd.DataFrame] | pd.DataFrame,
    purity_cutoff: float,
    total_counts_cutoff: int,
) -> pd.DataFrame:
    """Get the ZOTU count table for pure isolates in an isolate sequencing experiment.

    If df_purity is a list, we assume the dataframes have the same index, and a row
        is considered pure sa long as it is pure in ANY dataframes.

    Importantly, this function returns a binary count matrix whose index is from the
        index of `df_purity` and columns are from `top_1` column of `df_purity`.
        This is rather suboptimal. I need to think of an alternative way.
    """
    if isinstance(df_purity, list):
        df_purity = pd.concat(df_purity[["top_1", "top_1_val", "total_counts"]], axis=1)
    else:
        df_purity = df_purity[["top_1", "top_1_val", "total_counts"]]
    df_purity = df_purity.iloc[
        pd.concat(
            [
                (df_purity.iloc[:, i + 1] >= purity_cutoff)
                & (df_purity.iloc[:, i + 2] >= total_counts_cutoff)
                for i in range(0, df_purity.shape[1], 3)
            ],
            axis=1,
        )
        .any(axis=1)
        .tolist()
    ]
    top_1 = df_purity[["top_1"]].to_numpy()[
        np.arange(df_purity.shape[0]),
        df_purity[["top_1_val"]].to_numpy().argmax(axis=1),
    ]
    top_1 = pd.DataFrame(top_1, index=df_purity.index, columns=["#OTU ID"])
    return (
        top_1.reset_index(names="well_barcode")
        .assign(dummy=1)
        .pivot_table(
            index="well_barcode",
            columns="#OTU ID",
            values="dummy",
            fill_value=0,
        )
    ).astype(int)


def main():
# if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Utility functions for processing isolate metadata"
    )
    subparser = parser.add_subparsers(dest="subcommand")

    parser_isolate_meta = subparser.add_parser(
        "get_metadata",
        description="Summarize isolate metadata from the log files of isolate picking software",
    )
    parser_isolate_meta.add_argument(
        "-i",
        "--isolate_info_dir",
        help="Input directory containing the log file from isolate picking software",
        required=True,
        type=str,
        metavar="ISOLATE_METADATA_DIR",
    )
    parser_isolate_meta.add_argument(
        "-pm",
        "--plate_metadata_csv",
        help="Plate metadata CSV file",
        required=True,
        type=str,
        metavar="PLATE_METADATA_CSV",
    )
    parser_isolate_meta.add_argument(
        "-o",
        "--output_tsv",
        help="Output TSV file",
        required=True,
        type=str,
        metavar="OUTPUT_TSV",
    )

    parser_combine = subparser.add_parser(
        "combine_count",
        description="Combine zotu count tables",
    )
    parser_combine.add_argument(
        "-m",
        "--isolate_meta",
        help="Isolate metadata TSV file",
        required=True,
        type=str,
        metavar="ISOLATE_METADATA_TSV",
    )
    parser_combine.add_argument(
        "-cb",
        "--zotu_count_bacteria",
        help="Zotu count bacteria TSV file",
        default=None,
        type=str,
        metavar="ZOTU_COUNT_BACTERIA_TSV",
    )
    parser_combine.add_argument(
        "-cf",
        "--zotu_count_fungi",
        help="Zotu count fungi TSV file",
        default=None,
        type=str,
        metavar="ZOTU_COUNT_FUNGI_TSV",
    )
    parser_combine.add_argument(
        "-tb",
        "--taxonomy_bacteria",
        help="Taxonomy bacteria TSV file",
        default=None,
        type=str,
        metavar="TAXONOMY_BACTERIA_TSV",
    )
    parser_combine.add_argument(
        "-tf",
        "--taxonomy_fungi",
        help="Taxonomy fungi TSV file",
        default=None,
        type=str,
        metavar="TAXONOMY_FUNGI_TSV",
    )
    parser_combine.add_argument(
        "-pl",
        "--purity_levels",
        help="Taxonomic levels to calculate purity",
        default=[],
        type=str,
        nargs="+",
    )
    parser_combine.add_argument(
        "-cp",
        "--purity_cutoff",
        help="Purity cutoff",
        default=0.8,
        type=float,
    )
    parser_combine.add_argument(
        "-ct",
        "--total_counts_cutoff",
        help="Total counts cutoff",
        default=10,
        type=int,
    )
    parser_combine.add_argument(
        "-o",
        "--output_dir",
        help="Output directory",
        required=True,
        type=str,
        metavar="OUTPUT_DIR",
    )

    args = parser.parse_args()

    if args.subcommand == "get_metadata":
        isolate_info_dir = args.isolate_info_dir
        plate_metadata_csv = args.plate_metadata_csv
        output_tsv = args.output_tsv
        read_isolate_metadata_rich(isolate_info_dir, plate_metadata_csv).to_csv(
            output_tsv, sep="\t"
        )
    elif args.subcommand == "combine_count":
        (df_count_b, df_count_f), (df_taxon_b, df_taxon_f) = combine_count(
            args.isolate_meta,
            args.zotu_count_bacteria,
            args.zotu_count_fungi,
            args.taxonomy_bacteria,
            args.taxonomy_fungi,
        )
        os.makedirs(args.output_dir, exist_ok=True)
        if args.purity_levels:
            for l, (
                df_purity_bacteria,
                df_purity_fungi,
                df_purity,
            ) in get_isolate_purity(
                df_count_b,
                df_count_f,
                df_taxon_b,
                df_taxon_f,
                levels=args.purity_levels,
            ).items():
                if df_purity_bacteria is not None:
                    write_table(
                        df_purity_bacteria,
                        os.path.join(args.output_dir, f"purity_{l}_bacteria_only.tsv"),
                    )
                    otu_count_pure_b = (
                        quality_control(
                            df_purity_bacteria,
                            args.purity_cutoff,
                            args.total_counts_cutoff,
                        )
                        .reindex(df_count_b.index, fill_value=0)
                        .reindex(df_count_b.columns, axis=1, fill_value=0)
                    )
                    write_table(
                        otu_count_pure_b.transpose(),
                        os.path.join(args.output_dir, f"count_{l}_bacteria_only.biom"),
                    )
                if df_purity_fungi is not None:
                    write_table(
                        df_purity_fungi,
                        os.path.join(args.output_dir, f"purity_{l}_fungi_only.tsv"),
                    )
                    otu_count_pure_f = (
                        quality_control(
                            df_purity_fungi,
                            args.purity_cutoff,
                            args.total_counts_cutoff,
                        )
                        .reindex(df_count_f.index, fill_value=0)
                        .reindex(df_count_f.columns, axis=1, fill_value=0)
                    )
                    write_table(
                        otu_count_pure_f.transpose(),
                        os.path.join(args.output_dir, f"count_{l}_fungi_only.biom"),
                    )
                if df_purity is not None:
                    write_table(
                        df_purity, os.path.join(args.output_dir, f"purity_{l}.tsv")
                    )
                    otu_count_pure = (
                        quality_control(
                            df_purity, args.purity_cutoff, args.total_counts_cutoff
                        )
                        .reindex(df_count_b.index, fill_value=0)
                        .reindex(
                            df_count_b.columns.tolist() + df_count_f.columns.tolist(),
                            axis=1,
                            fill_value=0,
                        )
                        .transpose()
                    )
                    write_table(
                        otu_count_pure, os.path.join(args.output_dir, f"count_{l}.biom")
                    )
                    write_table(
                        otu_count_pure.loc[df_count_b.columns],
                        os.path.join(args.output_dir, f"count_{l}_bacteria.biom"),
                    )
                    write_table(
                        otu_count_pure.loc[df_count_f.columns],
                        os.path.join(args.output_dir, f"count_{l}_fungi.biom"),
                    )
        else:
            # simply concat them without any further modification, write data in the
            # same format as vsearch output
            df_zotu_count = (
                pd.concat([df_count_b, df_count_f], axis=1).fillna(0).transpose()
            )
            df_zotu_count.index.name = "#OTU ID"
            df_taxonomy = pd.concat([df_taxon_b, df_taxon_f], axis=0)
            write_table(df_zotu_count, os.path.join(args.output_dir, "count.biom"))
            df_taxonomy.to_csv(os.path.join(args.output_dir, "taxonomy.tsv"), sep="\t")
    else:
        raise ValueError("Invalid subcommand")
