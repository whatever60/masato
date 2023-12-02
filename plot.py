#!/usr/bin/env python3
import argparse
import os

import numpy as np
import pandas as pd
from skbio.diversity.alpha import chao1, shannon, simpson
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator


matplotlib.use("TkAgg")


def _get_color_palette(
    df: pd.DataFrame, group_key: str, column_of_interest: int, cmap: str = "viridis"
) -> dict:
    """Take a pandas dataframe, and return a dictionary mapping group to color. We
    expect the column of interest to be a continuous variable, and we will map each
    group to a color by its mean value of the column of interest.
    """
    mean = meta.groupby(group_key, sort=False)[column_of_interest].mean().reset_index()
    mean = mean.sort_values(by="richness")
    palette = sns.color_palette(cmap, n_colors=len(mean))
    color_map = {group: palette[i] for i, group in enumerate(mean[group_key])}
    return color_map


def _stacked_bar(
    dfs: list[pd.DataFrame],
    names: list[str],
    title: str,
    palette: list,
    axs: list[plt.Axes],
) -> None:
    for i, (ax, df, name) in enumerate(zip(axs, dfs, names)):
        # df = df.copy()
        # index = []
        # for i in df.index:
        #     if i.endswith("-Swab-combined"):
        #         index.append(i[1] + "-Bulk")
        #     elif i.endswith("-Scrape-R2A"):
        #         index.append(i[1] + "-Plate-R2A")
        #     elif i.endswith("-Scrape-TSA"):
        #         index.append(i[1] + "-Plate-TSA")
        #     else:
        #         raise ValueError(f"Unknown column: {i}")
        # df.index = index

        # df.index = [i for i in df.index]

        df.plot(kind="bar", stacked=True, color=palette, ax=ax)
        ax.set_xlabel("")
        ax.set_ylabel("Relative abundance")
        ax.set_title(name)
        # set xticklabels size
        # ax.tick_params(axis="x", labelsize=6)
        # remove legend unless the last plot
        if i != len(axs) - 1:
            ax.get_legend().remove()
        else:
            ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5))
            ax.get_legend().set_title(title)
        ax.set_ylabel(ax.yaxis.get_label().get_text(), fontsize=14)
        ax.set_xlabel(ax.xaxis.get_label().get_text(), fontsize=14)


def _heatmap(
    dfs: list[pd.DataFrame],
    names: list[str],
    title: str,
    cbar_label: str,
    axs: list[plt.Axes],
    cmap: str,
) -> None:
    for i, (ax, df, name) in enumerate(zip(axs, dfs, names)):
        # df = df.copy()
        # columns = []
        # for i in df.columns:
        #     if i.endswith("-Swab-combined"):
        #         columns.append(i[1] + "-Bulk")
        #     elif i.endswith("-Scrape-R2A"):
        #         columns.append(i[1] + "-Plate-R2A")
        #     elif i.endswith("-Scrape-TSA"):
        #         columns.append(i[1] + "-Plate-TSA")
        #     else:
        #         raise ValueError(f"Unknown column: {i}")
        # df.columns = columns

        if i == len(axs) - 1:
            cbar_ax = fig.add_axes(
                [
                    ax.get_position().x1 + 0.03,
                    axs[i - 1].get_position().y0,
                    0.03,
                    axs[i - 1].get_position().height,
                ]
            )
            sns.heatmap(
                df,
                ax=ax,
                cmap=cmap,
                lw=0.7,
                cbar=True,
                cbar_ax=cbar_ax,
                cbar_kws={"label": cbar_label},
                square=True,
            )
            cbar_ax.yaxis.label.set_size(12)
        sns.heatmap(
            df,
            ax=ax,
            cmap=cmap,
            lw=0.7,
            cbar=False,
            cbar_ax=None,
            square=True,
        )
        ax.set_xlabel("")
        ax.set_title(name)
        if i == 0:
            ax.set_ylabel(title, fontsize=16)
        # ax.tick_params(axis="x", labelsize=7)


def _barplot_with_whisker_strip(
    dfs: pd.DataFrame,
    names: list[str],
    title: str,
    group_key: str,
    column_of_interest: str,
    plot_strip: bool,
    axs: list[plt.Axes],
    ylog: bool = False,
) -> None:
    for i, (ax, df, name) in enumerate(zip(axs, dfs, names)):
        # df = df.copy()
        # index = []
        # for j in df[group_key]:
        #     if j.endswith("-Swab-combined"):
        #         index.append(j[1] + "-Bulk")
        #     elif j.endswith("-Scrape-R2A"):
        #         index.append(j[1] + "-Plate-R2A")
        #     elif j.endswith("-Scrape-TSA"):
        #         index.append(j[1] + "-Plate-TSA")
        #     else:
        #         raise ValueError(f"Unknown column: {j}")
        # df[group_key] = index

        group = df.groupby(group_key, sort=False)
        # Determine which groups have more than one sample
        groups_with_multiple_samples = group.filter(lambda x: len(x) > 1)[
            group_key
        ].unique()

        # assign different colors depending on if group.groups.keys() have "old", "reamplify", or "Scrape"
        colors = []
        for key in group.groups.keys():
            if key.endswith("-Bulk"):
                colors.append("tab:blue")
            elif key.endswith("-Plate-R2A"):
                colors.append("tab:orange")
            elif key.endswith("-Plate-TSA"):
                colors.append("tab:green")
            else:
                raise ValueError(f"Unknown group: {key}")

        # Plot bars with error bars for groups with multiple samples
        if len(groups_with_multiple_samples):
            ax.bar(
                range(len(group)),
                group[column_of_interest].mean(),
                yerr=group[column_of_interest].std()
                * [
                    1 if group in groups_with_multiple_samples else 0
                    for group in group.groups.keys()
                ],
                capsize=2,
                error_kw={"elinewidth": 1.5, "capthick": 1.5},
                width=0.6,
                # fill=False,
                color=colors,
            )

            # Remove edgecolor from bars to get rid of the faint gray line
            # for bar in bars:
            #     bar.set_edgecolor("none")

            # Use stripplot for individual points for groups with multiple samples
            if plot_strip:
                df_filtered = df[df[group_key].isin(groups_with_multiple_samples)]
                sns.stripplot(
                    x=group_key,
                    y=column_of_interest,
                    data=df_filtered,
                    color="gray",
                    order=group.groups.keys(),
                    ax=ax,
                    size=2,
                )
        else:
            ax.bar(
                range(len(group)),
                group[column_of_interest].mean(),
                width=0.5,
                # fill=False,
                color=colors,
            )
        ax.xaxis.set_major_locator(FixedLocator(range(len(group.groups.keys()))))
        ax.set_xticklabels(labels=group.groups.keys())
        ax.tick_params(axis="x", labelsize=8)

        ax.set_xlabel("")
        if i == 0:
            ax.set_ylabel(title, fontsize=14)
        else:
            ax.set_ylabel("")
        if ylog:
            ax.set_yscale("log")
        ax.set_title(name)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90)


def _calc_alpha_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Calculate alpha diversity metrics for each sample.
    Note that you should ignore chao1 if your input data is not integer count.
    """
    return pd.DataFrame(
        {
            "chao1": df.apply(chao1, axis=1),
            "richness": df.apply(lambda x: (x > 0).sum(), axis=1),
            "shannon": df.apply(shannon, axis=1),
            "simpson": df.apply(simpson, axis=1),
        },
        index=df.index,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command")

    parser_ab = subparsers.add_parser("abundance_group")
    parser_ab.add_argument(
        "plot_type",
        type=str,
        choices=["stacked_bar", "heatmap", "heatmap_log10", "all"],
    )
    parser_ab.add_argument("-f", "--fig_dir", type=str)
    parser_ab.add_argument("-m", "--metadata", type=str)
    parser_ab.add_argument("-i", "--rel_ab_dir", type=str)
    parser_ab.add_argument("-s", "--sample_group_key", type=str, default="sample_group")
    parser_ab.add_argument("-r", "--rep_group_key", type=str, default="rep_group")
    parser_ab.add_argument(
        "-t", "--tax_levels", nargs="+", default=["order", "family", "genus", "otu"]
    )

    parser_stats = subparsers.add_parser("stats_sample_count")
    parser_stats.add_argument("-f", "--fig_dir", type=str)
    parser_stats.add_argument("-i", "--abs_ab_dir", type=str)
    parser_stats.add_argument("-m", "--metadata", type=str)
    parser_stats.add_argument(
        "-s", "--sample_group_key", type=str, default="sample_group"
    )
    parser_stats.add_argument("-r", "--rep_group_key", type=str, default="rep_group")
    parser_stats.add_argument("-t", "--tax_levels", nargs="+", default=["otu"])
    parser_stats.add_argument("-p", "--plot_strip", action="store_true", default=False)

    parser_pcoa = subparsers.add_parser("dm")
    parser_pcoa.add_argument(
        "-i",
        "--input_ab",
        type=str,
        help="Input abundance matrix, in csv format, could be anything really, absolute "
        "or relative, sample-level or group-level. But probably sample x OTU count "
        "matrix is the most common one.",
    )
    parser_pcoa.add_argument("-m", "--metadata", type=str)
    parser_pcoa.add_argument("-f", "--fig_path", type=str)
    parser_pcoa.add_argument(
        "-s", "--sample_group_key", type=str, default="sample_group"
    )
    parser_pcoa.add_argument("-d", "--distance", type=str, choices=["bc", "euclid"])
    parser_pcoa.add_argument("-l", "--log10", action="store_true", default=False)
    parser_pcoa.add_argument(
        "-e", "--plot_ellipses", action="store_true", default=False
    )

    args = parser.parse_args()
    sample_group_key = args.sample_group_key
    metadata = args.metadata
    meta = pd.read_csv(metadata, sep=None, engine="python")

    if args.command in ["abundance_group", "stats_sample_count"]:
        fig_dir = args.fig_dir
        rep_group_key = args.rep_group_key
        tax_levels = args.tax_levels

        os.makedirs(fig_dir, exist_ok=True)
        names, groups = zip(*[i for i in meta.groupby(sample_group_key)])
        ratio = [i[rep_group_key].nunique() for i in groups]

        if args.command == "abundance_group":
            rel_ab_dir = args.rel_ab_dir
            plot_type = args.plot_type

            for level in tax_levels:
                res = pd.read_csv(f"{rel_ab_dir}/rel_ab_group_{level}.csv", index_col=0)
                res = res[sorted(res.columns)]
                res_group_list = [
                    res.loc[group[rep_group_key].unique()] for group in groups
                ]
                num_cols = len(res_group_list)
                width = 12
                wspace = 0.1

                if plot_type in ["stacked_bar", "all"]:
                    # stacked bar plot
                    fig_size = (width, 4)
                    custom_palette = (
                        sns.color_palette("tab20", 20)
                        + sns.color_palette("tab20b", 20)
                        + sns.color_palette("tab20c", 20)
                    )
                    fig, axs = plt.subplots(
                        1, num_cols, sharey=True, width_ratios=ratio, figsize=fig_size
                    )
                    _stacked_bar(
                        res_group_list,
                        names,
                        f"Taxonomy at {level if level != 'otu' else level.upper()} level",
                        custom_palette,
                        axs,
                    )
                    fig.subplots_adjust(wspace=wspace)
                    fig.savefig(
                        f"{fig_dir}/rel_ab_group_{level}_sb.png",
                        bbox_inches="tight",
                        dpi=300,
                    )
                    fig.savefig(
                        f"{fig_dir}/rel_ab_group_{level}_sb.pdf",
                        bbox_inches="tight",
                        dpi=300,
                    )
                if plot_type in ["heatmap", "heatmap_log10", "all"]:
                    size = (
                        res.shape[0] // (res.shape[0] / width),
                        res.shape[1] // (res.shape[0] / width),
                    )
                    if plot_type in ["heatmap_log10", "all"]:
                        fig, axs = plt.subplots(
                            1, num_cols, sharey=True, width_ratios=ratio, figsize=size
                        )
                        _heatmap(
                            [
                                np.log10(res_group + 1e-4).T
                                for res_group in res_group_list
                            ],
                            names,
                            title=f"Taxonomy at {level} level",
                            cbar_label="log10(relative abundance + 1e-4)",
                            axs=axs,
                            cmap="rocket_r",
                        )
                        fig.subplots_adjust(wspace=wspace)
                        fig.savefig(
                            f"{fig_dir}/rel_ab_group_{level}_hml.png",
                            bbox_inches="tight",
                            dpi=300,
                        )
                        fig.savefig(
                            f"{fig_dir}/rel_ab_group_{level}_hml.pdf",
                            bbox_inches="tight",
                            dpi=300,
                        )
                    if plot_type in ["heatmap", "all"]:
                        fig, axs = plt.subplots(
                            1,
                            len(groups),
                            sharey=True,
                            width_ratios=ratio,
                            figsize=size,
                        )
                        _heatmap(
                            [res_group.T for res_group in res_group_list],
                            names,
                            title=f"Taxonomy at {level} level",
                            cbar_label="relative abundance",
                            axs=axs,
                            cmap=None,
                        )
                        fig.subplots_adjust(wspace=wspace)
                        fig.savefig(
                            f"{fig_dir}/rel_ab_group_{level}_hm.png",
                            bbox_inches="tight",
                            dpi=300,
                        )
                        fig.savefig(
                            f"{fig_dir}/rel_ab_group_{level}_hm.pdf",
                            bbox_inches="tight",
                            dpi=300,
                        )
        elif args.command == "stats_sample_count":
            abs_ab_dir = args.abs_ab_dir
            plot_strip = args.plot_strip

            for level in tax_levels:
                res = pd.read_csv(f"{abs_ab_dir}/count_sample_{level}.csv", index_col=0)
                meta_l = pd.merge(
                    meta, _calc_alpha_metrics(res), left_on="sample", right_index=True
                )
                _, groups = zip(*[i for i in meta_l.groupby(sample_group_key)])
                title_prefix = (
                    f"{level.upper() if level == 'otu' else level.capitalize()}"
                )
                for column_of_interest, title in zip(
                    ["shannon", "simpson", "richness", "chao1", "read_count"],
                    [
                        f"{title_prefix} Shannon entropy",
                        f"{title_prefix} Simpson's index",
                        f"{title_prefix} Richness",
                        f"{title_prefix} Chao1 index",
                        "Library size",
                    ],
                ):
                    fig, axs = plt.subplots(
                        1, len(groups), sharey=True, width_ratios=ratio, figsize=(15, 3)
                    )
                    _barplot_with_whisker_strip(
                        groups,
                        names=names,
                        title=title,
                        group_key=rep_group_key,
                        column_of_interest=column_of_interest,
                        plot_strip=plot_strip,
                        axs=axs,
                        ylog=column_of_interest in ["read_count"],
                    )
                    fig.subplots_adjust(wspace=0.2)
                    fig.savefig(
                        f"{fig_dir}/{column_of_interest}_sample_{level}.png",
                        bbox_inches="tight",
                        dpi=300,
                    )
                    fig.savefig(
                        f"{fig_dir}/{column_of_interest}_sample_{level}.pdf",
                        bbox_inches="tight",
                        dpi=300,
                    )
    elif args.command == "dm":
        from plot_dm import plot_dm

        plot_dm(
            args.input_ab,
            args.metadata,
            args.fig_path,
            args.sample_group_key,
            args.distance,
            args.log10,
            args.plot_ellipses,
        )
    else:
        raise ValueError("Unknown command.")
