import argparse
import os

import numpy as np
import pandas as pd
from skbio.diversity.alpha import chao1, shannon, simpson
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator


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
        df.plot(kind="bar", stacked=True, color=palette, ax=ax)
        ax.set_xlabel("")
        ax.set_ylabel("Relative abundance")
        ax.set_title(name)
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


def _barplot_with_whisker_strip(
    dfs: pd.DataFrame,
    names: list[str],
    title: str,
    group_key: str,
    column_of_interest: str,
    plot_strip: bool,
    axs: list[plt.Axes],
) -> None:
    for i, (ax, df, name) in enumerate(zip(axs, dfs, names)):
        group = df.groupby(group_key, sort=False)
        # Determine which groups have more than one sample
        groups_with_multiple_samples = group.filter(lambda x: len(x) > 1)[
            group_key
        ].unique()

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
                capsize=3,
                error_kw={"elinewidth": 1.5, "capthick": 1.5},
                fill=False,
                width=0.6,
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
                fill=False,
                width=0.5,
            )
        ax.xaxis.set_major_locator(FixedLocator(range(len(group.groups.keys()))))
        ax.set_xticklabels(labels=group.groups.keys())

        ax.set_xlabel("")
        if i == 0:
            ax.set_ylabel(title, fontsize=14)
        else:
            ax.set_ylabel("")
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
        "plot_type", type=str, choices=["stacked_bar", "heatmap", "heatmap_log10"]
    )
    parser_ab.add_argument("-f", "--fig_dir", type=str)
    parser_ab.add_argument("-m", "--metadata", type=str)
    parser_ab.add_argument("-i", "--rel_ab_dir", type=str)
    parser_ab.add_argument("-s", "--sample_group_key", type=str, default="sample_group")
    parser_ab.add_argument("-r", "--rep_group_key", type=str, default="rep_group")
    parser_ab.add_argument(
        "-t", "--tax_levels", nargs="+", default=["order", "family", "genus", "otu"]
    )

    parser_stats = subparsers.add_parser("stats_sample")
    parser_stats.add_argument("-f", "--fig_dir", type=str)
    parser_stats.add_argument("-i", "--abs_ab_dir", type=str)
    parser_stats.add_argument("-m", "--metadata", type=str)
    parser_stats.add_argument(
        "-s", "--sample_group_key", type=str, default="sample_group"
    )
    parser_stats.add_argument("-r", "--rep_group_key", type=str, default="rep_group")
    parser_stats.add_argument("-t", "--tax_levels", nargs="+", default=["otu"])
    parser_stats.add_argument("-p", "--plot_strip", action="store_true", default=False)

    args = parser.parse_args()

    fig_dir = args.fig_dir
    metadata = args.metadata
    sample_group_key = args.sample_group_key
    rep_group_key = args.rep_group_key
    tax_levels = args.tax_levels

    os.makedirs(fig_dir, exist_ok=True)
    meta = pd.read_csv(metadata)
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

            if plot_type == "stacked_bar":
                # stacked bar plot
                custom_palette = (
                    sns.color_palette("tab20", 20)
                    + sns.color_palette("tab20b", 20)
                    + sns.color_palette("tab20c", 20)
                )
                fig, axs = plt.subplots(
                    1,
                    len(res_group_list),
                    sharey=True,
                    width_ratios=ratio,
                    figsize=(11, 4),
                )
                _stacked_bar(
                    res_group_list,
                    names,
                    f"Taxonomy at {level} level",
                    custom_palette,
                    axs,
                )
                fig.subplots_adjust(wspace=0.1)
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
            else:
                size = (
                    res.shape[0] // (res.shape[0] // 11),
                    res.shape[1] // (res.shape[0] // 11),
                )
                fig, axs = plt.subplots(
                    1, len(groups), sharey=True, width_ratios=ratio, figsize=size
                )
                if plot_type == "heatmap_log10":
                    _heatmap(
                        [np.log10(res_group + 1e-4).T for res_group in res_group_list],
                        names,
                        title=f"Taxonomy at {level} level",
                        cbar_label="log10(relative abundance + 1e-4)",
                        axs=axs,
                        cmap="coolwarm",
                    )
                    fig.subplots_adjust(wspace=0.3)
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
                elif plot_type == "heatmap":
                    _heatmap(
                        [res_group.T for res_group in res_group_list],
                        names,
                        title=f"Taxonomy at {level} level",
                        cbar_label="relative abundance",
                        axs=axs,
                        cmap=None,
                    )
                    fig.subplots_adjust(wspace=0.3)
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

    elif args.command == "stats_sample":
        abs_ab_dir = args.abs_ab_dir
        plot_strip = args.plot_strip

        for level in tax_levels:
            res = pd.read_csv(f"{abs_ab_dir}/abs_ab_sample_{level}.csv", index_col=0)
            meta_l = pd.merge(
                meta, _calc_alpha_metrics(res), left_on="sample", right_index=True
            )
            _, groups = zip(*[i for i in meta_l.groupby(sample_group_key)])
            for column_of_interest, title in zip(
                ["shannon", "simpson", "richness", "chao1"],
                ["Shannon entropy", "Simpson's index", "OTU richness", "Chao1 index"],
            ):
                fig, axs = plt.subplots(
                    1, len(groups), sharey=True, width_ratios=ratio, figsize=(15, 3)
                )
                _barplot_with_whisker_strip(
                    groups,
                    names=names,
                    title=f"{level.upper() if title == 'otu' else level.capitalize()} {title}",
                    group_key=rep_group_key,
                    column_of_interest=column_of_interest,
                    plot_strip=plot_strip,
                    axs=axs,
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
    else:
        raise ValueError("Unknown command.")
