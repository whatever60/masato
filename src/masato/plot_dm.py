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

from masato.get_abundance import get_otu_count, _agg_along_axis


def _load(ab_path: str, metadata_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_table(metadata_path, index_col=0)
    df = pd.read_table(ab_path, index_col=0).T

    # drop_samples = df.index.str.contains("-V-") | df.index.str.contains("-SNS-")
    # drop_samples = df.index.str.contains("-5N-")
    # df = df.loc[~drop_samples]
    # meta = meta.loc[~drop_samples]

    nonzero = df.sum(axis=1).astype(bool)
    df = df.loc[nonzero]
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
    fig_path: str | None = None,
    distance: str = "braycurtis",
    title: str | None = None,
    hue: str | None = None,
    hue_order: list[str] | None = None,
    hue_dict: str | dict | None = None,
    style: str | None = None,
    style_order: list[str] | None = None,
    style_dict: dict[str] | bool = True,
    s: int = 50,
    annotate_dots: bool = False,
    plot_ellipses: bool = False,
    plot_pc3: bool = True,
    fig: matplotlib.figure.Figure | None = None,
    axs: list[matplotlib.axes.Axes] | matplotlib.axes.Axes | None = None,
) -> None:
    n_components = min(10, *df_otu_count.shape)
    n_pcs = min(3, n_components)
    if n_components == 1:
        raise ValueError(
            f"Minimum dimensionality of input data is 2, getting {n_components}."
        )

    if distance == "braycurtis":
        bd = beta_diversity("braycurtis", df_otu_count)
        dist_mtx = np.nan_to_num(bd.data, 0)
        pc_obj = pcoa(dist_mtx, number_of_dimensions=n_components)
        pc = pc_obj.samples.copy()
        pc.index = df_otu_count.index
        variance = pc_obj.proportion_explained.to_numpy()
    elif distance == "euclid":
        pc_obj = PCA(n_components=n_components).fit(df_otu_count)
        pc = pc_obj.transform(df_otu_count)
        pc = pd.DataFrame(
            pc,
            index=df_otu_count.index,
            columns=[f"PC{i}" for i in range(1, pc.shape[1] + 1)],
        )
        variance = pc_obj.explained_variance_ratio_
    else:
        raise NotImplementedError(f"Unsupported distance metric: {distance}")

    if n_pcs == 2:
        pc["PC3"] = 0
        variance = np.append(variance, 0)

    marker_size = 8
    pc = pd.concat([pc, df_meta], axis=1)

    if hue is None:
        pc["_hue"] = "all"
        hue = "_hue"
    if hue_order is None:
        hue_order = sorted(pc[hue].unique().tolist())

    if style is None:
        style_order = None
    elif style_order is None:
        style_order = sorted(pc[style].unique().tolist())

    axis_label_fs = 14
    title_fs = 16

    # prepare figure and axes
    if fig is None or axs is None:
        if plot_pc3:
            fig, axs = plt.subplots(1, 2, figsize=(8, 3.5))
        else:
            fig, ax = plt.subplots(1, 1, figsize=(5.5, 4))
            axs = [ax]
    else:
        expected_axs = 2 if plot_pc3 else 1
        if isinstance(axs, matplotlib.axes.Axes):
            axs = [axs]
        if len(axs) != expected_axs:
            raise ValueError(f"Expected {expected_axs} axes, got {len(axs)}.")

    axs: list[matplotlib.axes.Axes]
    sns.scatterplot(
        data=pc,
        x="PC1",
        y="PC2",
        hue=hue,
        hue_order=hue_order,
        palette=hue_dict,
        style=style,
        markers=style_dict,
        style_order=style_order,
        linewidth=0,
        s=s,
        legend=not plot_pc3,  # Show legend here only if it's the single plot
        ax=axs[0],
    )
    axs[0].set_xlabel(f"PC1 ({variance[0] * 100:.2f}%)", fontsize=axis_label_fs)
    axs[0].set_ylabel(f"PC2 ({variance[1] * 100:.2f}%)", fontsize=axis_label_fs)

    if plot_pc3:
        sns.scatterplot(
            data=pc,
            x="PC2",
            y="PC3",
            hue=hue,
            hue_order=hue_order,
            palette=hue_dict,
            style=style,
            markers=style_dict,
            style_order=style_order,
            linewidth=0,
            s=s,
            legend=True,  # Legend is on this plot if it exists
            ax=axs[1],
        )
        axs[1].set_xlabel(f"PC2 ({variance[1] * 100:.2f}%)", fontsize=axis_label_fs)
        axs[1].set_ylabel(f"PC3 ({variance[2] * 100:.2f}%)", fontsize=axis_label_fs)

    # set axis to be square
    for ax in axs:
        xleft, xright = ax.get_xlim()
        ybottom, ytop = ax.get_ylim()
        ax.set_aspect(abs((xright - xleft) / (ybottom - ytop)), adjustable="box")
    legend_ax = axs[-1]  # Target the last axis for legend modifications
    if hue != "_hue":
        handles, labels = legend_ax.get_legend_handles_labels()
        hues = pc[hue].unique().tolist()
        for h, label in zip(handles, labels):
            if label in hues:
                h.set_marker("o")
            h.set_markersize(marker_size)
        legend_ax.legend(handles, labels)
    legend_ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0.0)

    if annotate_dots:
        for i, txt in enumerate(pc.index):
            axs[0].annotate(
                txt, (pc["PC1"].iloc[i], pc["PC2"].iloc[i]), fontsize=6, alpha=0.4
            )
            if plot_pc3:
                axs[1].annotate(
                    txt, (pc["PC2"].iloc[i], pc["PC3"].iloc[i]), fontsize=6, alpha=0.4
                )

    if plot_ellipses:
        pc_dims_list, ax_list = (
            ([(1, 2), (2, 3)], axs) if plot_pc3 else ([(1, 2)], [axs[0]])
        )
        n_hues = len(hue_order)
        palette = (
            sns.color_palette(hue_dict, n_hues)
            if isinstance(hue_dict, str)
            else (
                list(hue_dict.values())
                if isinstance(hue_dict, dict)
                else sns.color_palette(n_colors=n_hues)
            )
        )
        hue_to_color = dict(zip(hue_order, palette))

        for pc_dims, ax in zip(pc_dims_list, ax_list):
            for group in hue_order:
                group_data = pc.loc[pc[hue] == group]
                x = group_data[f"PC{pc_dims[0]}"]
                y = group_data[f"PC{pc_dims[1]}"]
                if len(x) < 2:
                    continue

                mean_x, mean_y = x.mean(), y.mean()
                cov = np.cov(x, y)
                lambda_, v = eigsorted(cov)
                theta = np.degrees(np.arctan2(*v[:, 0][::-1]))
                w, h = 4 * np.sqrt(lambda_)  # 2 * std_dev
                ell = Ellipse(
                    xy=(mean_x, mean_y), width=w, height=h, angle=theta, alpha=0.2
                )
                ell.set_facecolor(hue_to_color[group])
                ell.set_edgecolor("grey")
                ax.add_artist(ell)

    if title:
        fig.suptitle(title, fontsize=title_fs)

    fig.tight_layout()
    # Adjust layout after tight_layout to prevent title overlap
    if title:
        fig.subplots_adjust(top=0.88)

    if fig_path:
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        fig.savefig(os.path.splitext(fig_path)[0] + ".svg", bbox_inches="tight")

    return fig, axs


def main():
    # if __name__ == "__main__":
    import argparse

    # configure matplotlib PDF saving to use text instead of vector graphics
    plt.rcParams["pdf.fonttype"] = 42
    plt.rcParams["svg.fonttype"] = "none"
    # matplotlib.use("TkAgg")

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
        "-s",
        "--sample_group_key",
        default="sample_group",
        help="Column name in the metadata table to group samples",
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
        "-a",
        "--annotate_dots",
        help="Annotate dots with sample names",
        action="store_true",
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
    sample_group_key = args.sample_group_key
    rep_group_key = args.rep_group_key
    spikein_taxa_key = args.spikein_taxa_key
    sample_weight_key = args.sample_weight_key
    log10 = args.log10
    distance = args.distance
    hue = args.hue
    style = args.style
    annotate_dots = args.annotate_dots
    ellipses = args.ellipses

    df_otu_count, df_meta, df_tax = get_otu_count(
        otu_count_tsv,
        metadata,
        otu_taxonomy_tsv,
        sample_weight_key=sample_weight_key,
        spikein_taxa_key=spikein_taxa_key,
    )

    # [PCA|PCoA (with Bray-Curtis distance)] on[log10 | ]<[tax_level> [relative abundance|absolute abundance|raw count]
    # use keyword format
    title_template = "{distance} on {log10}{tax_level} {transform}"
    title_kw = {}
    if transform == "relative":
        df_otu = df_otu_count.div(df_otu_count.sum(axis=1), axis=0)
        title_kw["transform"] = "relative abundance"
    elif transform == "absolute":
        df_otu = df_otu_count.div(df_meta.norm_factor, axis=0)
        title_kw["transform"] = "absolute abundance"
    elif transform == "count":
        df_otu = df_otu_count.copy()
        title_kw["transform"] = "raw count"
    else:
        raise ValueError(f"Unsupported transform: {transform}")

    if sample_group_key is not None and rep_group_key is not None:
        group = df_meta[sample_group_key] + "-" + df_meta[rep_group_key]
        df_otu = _agg_along_axis(df_otu, group, axis=0)
        df_meta = df_meta.groupby(group).first()
    if log10:
        title_kw["log10"] = "log10 "
        df_otu = np.log10(df_otu + 1e-8)
    else:
        title_kw["log10"] = ""

    if distance == "braycurtis":
        title_kw["distance"] = "PCoA (Bray-Curtis)"
    elif distance == "euclid":
        title_kw["distance"] = "PCA"
    else:
        raise ValueError(
            f"Unsupported distance metric: {distance}, select from 'braycurtis' or 'euclid'."
        )

    os.makedirs(fig_dir, exist_ok=True)
    for level in tax_levels:
        df_otu_tax = _agg_along_axis(df_otu, df_tax[level], axis=1)
        fig_path = os.path.join(
            fig_dir,
            f"{distance}_{'log10_' if log10 else ''}{transform}_{'g' if rep_group_key is not None else 's'}_{level}.png",
        )
        if level == "otu":
            level = "ZOTU"
        else:
            level = level.capitalize()
        plot_dm(
            df_otu_tax,
            df_meta,
            fig_path,
            distance=distance,
            title=title_template.format(**title_kw, tax_level=level),
            hue=hue,
            style=style,
            annotate_dots=annotate_dots,
            plot_ellipses=ellipses,
        )
