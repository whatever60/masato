#!/usr/bin/env python
import argparse
import os
import glob

import numpy as np
import pandas as pd

from get_abundance import get_otu_count, _agg_along_axis, _taxa_qc


def _read_isolate_metadata(isolate_metadata_dir: str) -> pd.DataFrame:
    isolate_metadata_list = []
    for path in glob.glob(
        os.path.join(isolate_metadata_dir, "Destination * - * *.csv")
    ):
        dest_plate_barcode = os.path.basename(path).split(" ")[1]
        df = pd.read_csv(path, skiprows=2)
        df["dest_plate"] = dest_plate_barcode
        isolate_metadata_list.append(df)
    return pd.concat(isolate_metadata_list, ignore_index=True)


def read_isolate_metadata_rich(
    isolate_metadata_dir: str, plate_metadata_path: str
) -> pd.DataFrame:
    plate_metadata = pd.read_csv(plate_metadata_path)
    if not plate_metadata.barcode.is_unique:
        raise ValueError("Plate metadata barcode column is not unique")
    isolate_metadata = _read_isolate_metadata(isolate_metadata_dir)
    isolate_metadata.columns = ["picking_coord", "src_plate", "dest_well", "dest_plate"]
    isolate_metadata = isolate_metadata.merge(
        plate_metadata[
            ["barcode", "medium_type", "sample_type", "sample_group"]
        ].rename({"barcode": "src_plate"}, axis=1),
        on="src_plate",
        how="left",
    )
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
    return isolate_metadata.set_index("sample").sort_index()


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


def get_isolate_purity(
    df_isolate_meta_tsv: str,
    df_zotu_count_bacteria_tsv: str = None,
    df_zotu_count_fungi_tsv: str = None,
    df_taxonomy_bacteria_tsv: str = None,
    df_taxonomy_fungi_tsv: str = None,
    level: str = "otu",
):
    """Get purity table for isolates.

    If bacteria zotu count table is given, get its purity.
    If fungi zotu count table is given, get its purity.
    If both are given, get the purity for bacteria and fungi individually and collectively.

    If taxonomy tables are given, attach the taxonomy information to the purity table.
    """
    both = 0
    ret = []
    if df_zotu_count_bacteria_tsv is not None:
        df_zotu_cb, _, df_taxonomy_b = get_otu_count(
            df_zotu_count_bacteria_tsv,
            df_isolate_meta_tsv,
            df_taxonomy_bacteria_tsv,
            warning=False,
        )
        df_rab_b = _taxa_qc(
            _agg_along_axis(
                df_zotu_cb.div(df_zotu_cb.sum(axis=1), axis=0).fillna(0),
                df_taxonomy_b[level],
                axis=1,
            ),
            rel_ab_thres=0,
            keep_rare=True,
            keep_unknown=False,
        )
        df_purity_bacteria = get_top_cols(df_rab_b)
        df_purity_bacteria["total_counts"] = df_zotu_cb.sum(axis=1)
        both += 1
        ret.append(df_purity_bacteria)
    else:
        ret.append(None)
    if df_zotu_count_fungi_tsv is not None:
        df_zotu_cf, _, df_taxonomy_f = get_otu_count(
            df_zotu_count_fungi_tsv,
            df_isolate_meta_tsv,
            df_taxonomy_fungi_tsv,
            warning=False,
        )
        df_rab_f = _taxa_qc(
            _agg_along_axis(
                df_zotu_cf.div(df_zotu_cf.sum(axis=1), axis=0).fillna(0),
                df_taxonomy_f[level],
                axis=1,
            ),
            rel_ab_thres=0,
            keep_rare=True,
            keep_unknown=False,
        )
        df_purity_fungi = get_top_cols(df_rab_f)
        df_purity_fungi["total_counts"] = df_zotu_cf.sum(axis=1)
        both += 1
        ret.append(df_purity_fungi)
    else:
        ret.append(None)
    if both == 2:
        df_zotu_count = pd.concat([df_zotu_cb, df_zotu_cf], axis=1).fillna(0)
        df_taxonomy = pd.concat([df_taxonomy_b, df_taxonomy_f], axis=0)
        df_rab = _taxa_qc(
            _agg_along_axis(
                df_zotu_count.div(df_zotu_count.sum(axis=1), axis=0).fillna(0),
                df_taxonomy[level],
                axis=1,
            ),
            rel_ab_thres=0,
            keep_rare=True,
            keep_unknown=False,
        )
        df_purity = get_top_cols(df_rab)
        df_purity["total_counts"] = df_zotu_count.sum(axis=1)
        ret.append(df_purity)
    else:
        ret.append(None)
    [df.sort_index(inplace=True) for df in ret if df is not None]
    return ret


def combine_count(
    df_isolate_meta_tsv: str,
    df_zotu_count_bacteria_tsv: str,
    df_zotu_count_fungi_tsv: str,
    df_taxonomy_bacteria_tsv: str = None,
    df_taxonomy_fungi_tsv: str = None,
) -> pd.DataFrame:
    df_zotu_cb, _, df_taxonomy_bacteria = get_otu_count(
        df_zotu_count_bacteria_tsv,
        df_isolate_meta_tsv,
        df_taxonomy_bacteria_tsv,
        warning=False,
    )
    df_zotu_cf, _, df_taxonomy_fungi = get_otu_count(
        df_zotu_count_fungi_tsv,
        df_isolate_meta_tsv,
        df_taxonomy_fungi_tsv,
        warning=False,
    )
    return (
        pd.concat([df_zotu_cb, df_zotu_cf], axis=1).fillna(0),
        pd.concat([df_taxonomy_bacteria, df_taxonomy_fungi]),
    )


if __name__ == "__main__":
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

    parser_purity = subparser.add_parser(
        "get_purity",
        description="Get purity table for isolates",
    )
    parser_purity.add_argument(
        "-m",
        "--isolate_meta",
        help="Isolate metadata TSV file",
        required=True,
        type=str,
        metavar="ISOLATE_METADATA_TSV",
    )
    parser_purity.add_argument(
        "-cb",
        "--zotu_count_bacteria",
        help="Zotu count bacteria TSV file",
        default=None,
        type=str,
        metavar="ZOTU_COUNT_BACTERIA_TSV",
    )
    parser_purity.add_argument(
        "-cf",
        "--zotu_count_fungi",
        help="Zotu count fungi TSV file",
        default=None,
        type=str,
        metavar="ZOTU_COUNT_FUNGI_TSV",
    )
    parser_purity.add_argument(
        "-tb",
        "--taxonomy_bacteria",
        help="Taxonomy bacteria TSV file",
        default=None,
        type=str,
        metavar="TAXONOMY_BACTERIA_TSV",
    )
    parser_purity.add_argument(
        "-tf",
        "--taxonomy_fungi",
        help="Taxonomy fungi TSV file",
        default=None,
        type=str,
        metavar="TAXONOMY_FUNGI_TSV",
    )
    parser_purity.add_argument(
        "-l",
        "--level",
        help="Taxonomic level",
        default="otu",
        type=str,
        metavar="LEVEL",
    )
    parser_purity.add_argument(
        "-o",
        "--output_dir",
        help="Output directory",
        required=True,
        type=str,
        metavar="OUTPUT_DIR",
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
        required=True,
        type=str,
        metavar="ZOTU_COUNT_BACTERIA_TSV",
    )
    parser_combine.add_argument(
        "-cf",
        "--zotu_count_fungi",
        help="Zotu count fungi TSV file",
        required=True,
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
        "-oc",
        "--output_count_tsv",
        help="Output count TSV file",
        required=True,
        type=str,
        metavar="OUTPUT_COUNT_TSV",
    )
    parser_combine.add_argument(
        "-ot",
        "--output_taxonomy_tsv",
        help="Output taxonomy TSV file",
        required=False,
        type=str,
        metavar="OUTPUT_TAXONOMY_TSV",
    )

    args = parser.parse_args()

    if args.subcommand == "get_metadata":
        isolate_info_dir = args.isolate_info_dir
        plate_metadata_csv = args.plate_metadata_csv
        output_tsv = args.output_tsv
        read_isolate_metadata_rich(isolate_info_dir, plate_metadata_csv).to_csv(
            output_tsv, sep="\t"
        )
    elif args.subcommand == "get_purity":
        level = args.level
        df_purity_bacteria, df_purity_fungi, df_purity = get_isolate_purity(
            args.isolate_meta,
            args.zotu_count_bacteria,
            args.zotu_count_fungi,
            args.taxonomy_bacteria,
            args.taxonomy_fungi,
            level,
        )
        os.makedirs(args.output_dir, exist_ok=True)
        if df_purity_bacteria is not None:
            df_purity_bacteria.to_csv(
                os.path.join(args.output_dir, f"purity_{level}_bacteria.tsv"), sep="\t"
            )
        if df_purity_fungi is not None:
            df_purity_fungi.to_csv(
                os.path.join(args.output_dir, f"purity_{level}_fungi.tsv"), sep="\t"
            )
        if df_purity is not None:
            df_purity.to_csv(
                os.path.join(args.output_dir, f"purity_{level}.tsv"), sep="\t"
            )
    elif args.subcommand == "combine_count":
        combined_count, combined_tax = combine_count(
            args.isolate_meta,
            args.zotu_count_bacteria,
            args.zotu_count_fungi,
            args.taxonomy_bacteria,
            args.taxonomy_fungi,
        )
        combined_count.to_csv(args.output_count_tsv, sep="\t")
        combined_tax.to_csv(args.output_taxonomy_tsv, sep="\t")
    else:
        raise ValueError("Invalid subcommand")
