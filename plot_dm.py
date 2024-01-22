#!/usr/bin/env python3
"""
Figures that have been tested:
- Euclidean distance on log10 relative abundance at sample x genus (or others) level
- Bray-Curtis distance on count at sample x ZOTU (or others) level
"""

import os

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns

from get_abundance import get_otu_count, _agg_along_axis

# configure matplotlib PDF saving to use text instead of vector graphics
plt.rcParams["pdf.fonttype"] = 42
matplotlib.use("TkAgg")


def _load(ab_path: str, metadata_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_table(metadata_path, index_col=0)
    df = pd.read_table(ab_path, index_col=0).T

    # drop_samples = df.index.str.contains("-V-") | df.index.str.contains("-SNS-")
    # drop_samples = df.index.str.contains("-5N-")
    # df = df.loc[~drop_samples]
    # meta = meta.loc[~drop_samples]

    nonzero = df.sum(axis=1).astype(bool)
    df = df.loc[nonzero]
    # import pdb; pdb.set_trace()
    meta = meta.loc[nonzero]

    # source = []
    # for i in meta.index:
    #     if "-S-" in i:
    #         if "-OGAmp2-" in i:
    #             source.append("Swab-new")
    #         elif "-OG-" in i:
    #             source.append("Swab-old")
    #         elif "-OGcomb-" in i:
    #             source.append("Swab-comb")
    #         else:
    #             raise ValueError(i)
    #     if "-V-" in i:
    #         source.append("Vortex")
    #     if "-SNS-" in i:
    #         source.append("Swab no spike-in")
    #     if "-R-" in i:
    #         source.append("R2A")
    #     if "-T-" in i:
    #         source.append("TSA")
    # meta["source"] = source

    rep_number = []
    for i in meta.index:
        # if "-S-" in i:
        if i.startswith("B"):
            rep_number.append("AFRL Boneyard")
        elif i.startswith("C"):
            rep_number.append("Cell concentrate")
        else:
            raise ValueError(i)
    meta["Source"] = rep_number
    # meta["rep"] = "every_replication"

    meta["source"] = "every_source"
    return df, meta


def eigsorted(cov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


def plot_dm(
    df_otu_count: pd.DataFrame,
    df_meta: pd.DataFrame,
    fig_path: str,
    distance: str,
    hue: str = None,
    style: str = None,
    plot_ellipses: bool = False,
) -> None:
    if distance == "braycurtis":
        pc_obj = pcoa(
            beta_diversity("braycurtis", df_otu_count), number_of_dimensions=10
        )
        pc = pc_obj.samples.copy()
        pc.index = df_otu_count.index
        variance = pc_obj.proportion_explained.to_numpy()
    elif distance == "euclid":
        pc_obj = PCA(n_components=min(10, df_otu_count.shape[1])).fit(df_otu_count)
        pc = pc_obj.transform(df_otu_count)
        pc = pd.DataFrame(
            pc,
            index=df_otu_count.index,
            columns=[f"PC{i}" for i in range(1, pc.shape[1] + 1)],
        )
        variance = pc_obj.explained_variance_ratio_
    else:
        raise ValueError(
            f"Unsupported distance metric: {distance}, select from 'braycurtis' or 'euclid'."
        )

    marker_size = 8
    pc = pd.concat([pc, df_meta], axis=1)

    # markers = {"Bulk": "X", "Plate-R2A": "o", "Plate-TSA": "s"}
    if hue is None:
        pc["_hue"] = "all"
        hue = "_hue"
    if style is None:
        style_order = None
    else:
        style_order = sorted(pc[style].unique().tolist())

    fig, axs = plt.subplots(1, 2, figsize=(8, 3))
    sns.scatterplot(
        data=pc,
        x="PC1",
        y="PC2",
        hue=hue,
        style=style,
        linewidth=0,
        legend=False,
        ax=axs[0],
    )
    axs[0].set_xlabel(f"PC1 ({variance[0] * 100:.2f}%)")
    axs[0].set_ylabel(f"PC2 ({variance[1] * 100:.2f}%)")

    sns.scatterplot(
        data=pc,
        x="PC2",
        y="PC3",
        hue=hue,
        style=style,
        # markers=markers,
        style_order=style_order,
        linewidth=0,
        ax=axs[1],
    )
    axs[1].set_xlabel(f"PC2 ({variance[1] * 100:.2f}%)")
    axs[1].set_ylabel(f"PC3 ({variance[2] * 100:.2f}%)")

    if hue != "_hue":
        # set legend markers under `hue` to o
        hues = pc[hue].unique().tolist()
        handles, labels = axs[1].get_legend_handles_labels()
        for h, l in zip(handles, labels):
            if l in hues:
                # set marker size to s and marker to o
                h.set_marker("o")
            h.set_markersize(marker_size)
        axs[1].legend(handles, labels)

    axs[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

    if plot_ellipses:
        for pc_dims, ax in zip([(1, 2), (2, 3)], axs):
            n = len(pc[hue].unique())
            row_colors = sns.color_palette("tab10", n)
            for group in pc[hue].unique():
                x = pc.loc[pc[hue] == group, f"PC{pc_dims[0]}"]
                y = pc.loc[pc[hue] == group, f"PC{pc_dims[1]}"]
                mean_x, mean_y = x.mean(), y.mean()
                cov = np.cov(x, y)
                lambda_, v = eigsorted(cov)
                theta = np.degrees(np.arctan2(*v[:, 0][::-1]))
                w, h = 2 * 2 * np.sqrt(lambda_)
                ell = Ellipse(
                    xy=(mean_x, mean_y),
                    width=w,
                    height=h,
                    angle=theta,
                    # color=row_colors[pc["sample_group"].unique().tolist().index(group)],
                    alpha=0.2,
                )
                ell.set_facecolor(row_colors[pc[hue].unique().tolist().index(group)])
                ell.set_edgecolor("grey")
                ax.add_artist(ell)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    # also save a pdf file
    fig.savefig(os.path.splitext(fig_path)[0] + ".pdf", bbox_inches="tight")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--otu_count_tsv",
        help="Path to the OTU count table",
        required=True,
    )
    parser.add_argument(
        "-m",
        "--metadata",
        help="Path to the metadata table",
        required=True,
    )
    parser.add_argument("-t", "--otu_taxonomy_tsv", type=str, required=True)
    parser.add_argument(
        "-f",
        "--fig_dir",
        help="Path to the output figure",
        required=True,
    )
    parser.add_argument(
        "-l",
        "--tax_levels",
        help="Taxonomic levels to plot",
        nargs="+",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--rep_group_key",
        default="rep_group",
        help="Column name in the metadata table to group replicates. If null, plot is "
        "generated without aggregating replicates.",
    )
    parser.add_argument("-sp", "--spikein_taxa_key", type=str, default="spike_in")
    parser.add_argument("-w", "--sample_weight_key", default="sample_weight", type=str)
    # parser.add_argument(
    #     "--relative",
    #     help="Plot relative abundance",
    #     action="store_true",
    # )
    # parser.add_argument(
    #     "--absolute",
    #     help="Plot absolute abundance, exclusive with --relative",
    #     action="store_true",
    # )
    parser.add_argument(
        "-tr",
        "--transform",
        choices=["count", "relative", "absolute"],
        default="count",
        help="Transform to apply to the OTU count table",
    )
    parser.add_argument(
        "-lg",
        "--log10",
        help="Log10 transform",
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--distance",
        help="Distance metric",
        required=True,
        choices=["braycurtis", "euclid"],
    )
    parser.add_argument(
        "--hue",
        help="Column name in the metadata table to color samples",
        default=None,
    )
    parser.add_argument(
        "--style",
        help="Column name in the metadata table to style samples",
        default=None,
    )

    parser.add_argument(
        "-e",
        "--ellipses",
        help="Plot ellipses",
        action="store_true",
    )
    args = parser.parse_args()

    otu_count_tsv = args.otu_count_tsv
    metadata = args.metadata
    otu_taxonomy_tsv = args.otu_taxonomy_tsv
    fig_dir = args.fig_dir
    tax_levels = args.tax_levels
    transform = args.transform
    rep_group_key = args.rep_group_key
    spikein_taxa_key = args.spikein_taxa_key
    sample_weight_key = args.sample_weight_key
    log10 = args.log10
    distance = args.distance
    hue = args.hue
    style = args.style
    ellipses = args.ellipses

    df_otu_count, df_meta, df_tax = get_otu_count(
        otu_count_tsv,
        metadata,
        otu_taxonomy_tsv,
        sample_weight_key=sample_weight_key,
        spikein_taxa_key=spikein_taxa_key,
    )
    if transform == "relative":
        df_otu = df_otu_count.div(df_otu_count.sum(axis=1), axis=0)
    elif transform == "absolute":
        df_otu = df_otu_count.div(df_meta.norm_factor, axis=0)
    elif transform == "count":
        df_otu = df_otu_count.copy()
    else:
        raise ValueError(f"Unsupported transform: {transform}")

    if rep_group_key is not None:
        df_otu = _agg_along_axis(df_otu, df_meta[rep_group_key], axis=0)
        df_meta = df_meta.groupby(rep_group_key).first()
    if log10:
        df_otu = np.log10(df_otu + 1e-8)

    os.makedirs(fig_dir, exist_ok=True)
    for level in tax_levels:
        df_otu_tax = _agg_along_axis(df_otu, df_tax[level], axis=1)
        fig_path = os.path.join(
            fig_dir,
            f"{distance}_{'log10_' if log10 else ''}{transform}_{'g' if rep_group_key is not None else 's'}_{level}.png",
        )
        plot_dm(
            df_otu_tax,
            df_meta,
            fig_path,
            distance=distance,
            hue=hue,
            style=style,
            plot_ellipses=ellipses,
        )
