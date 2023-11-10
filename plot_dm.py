import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import seaborn as sns


def _load(ab_path: str, metadata_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = pd.read_csv(metadata_path, index_col=0)
    df = pd.read_csv(ab_path, index_col=0)

    # drop_samples = df.index.str.contains("-V-") | df.index.str.contains("-SNS-")
    drop_samples = df.index.str.contains("-5N-")
    df = df.loc[~drop_samples]
    meta = meta.loc[~drop_samples]
    
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

    # print(np.unique([i.split("-")[2] for i in meta.index]))
    rep1, rep2 = np.unique([i.split("-")[2] for i in meta.index]).reshape(-1, 2).T
    # print(rep1, rep2)
    rep_number = []
    for i in meta.index:
        if any(f"-{r}-" in i for r in rep1):
            rep_number.append("rep 1")
        elif any(f"-{r}-" in i for r in rep2):
            rep_number.append("rep 2")
        else:
            raise ValueError(i)
    meta["rep"] = rep_number
    # meta["rep"] = "every_replication"

    meta["source"] = "every_source"
    return df, meta


def eigsorted(cov: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    return vals[order], vecs[:, order]


def plot_dm(
    ab_path: str,
    meta_path: str,
    fig_path: str,
    sample_group_key: str,
    distance: str,
    log10: False,
    plot_ellipses: bool = False,
) -> None:
    ab, meta = _load(ab_path, meta_path)
    if log10:
        ab = np.log10(ab + 1e-8)
    if distance == "bc":
        pc_obj = pcoa(beta_diversity("braycurtis", ab), number_of_dimensions=10)
        pc = pc_obj.samples.copy()
        pc.index = ab.index
        variance = pc_obj.proportion_explained.to_numpy()
    elif distance == "euclid":
        pc_obj = PCA(n_components=10).fit(ab)
        pc = pc_obj.transform(ab)
        pc = pd.DataFrame(
            pc,
            index=ab.index,
            columns=[f"PC{i}" for i in range(1, pc.shape[1] + 1)],
        )
        variance = pc_obj.explained_variance_ratio_
    else:
        raise ValueError(
            f"Unsupported distance metric: {distance}, select from 'bc' or 'euclid'"
        )
    pc[["Sample group", "Source", "Replication number"]] = meta[
        [sample_group_key, "source", "rep"]
    ]

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    sns.scatterplot(
        data=pc,
        x="PC1",
        y="PC2",
        hue="Sample group",
        style="Replication number",
        ax=axs[0],
        legend=False,
    )
    axs[0].set_xlabel(f"PC1 ({variance[0] * 100:.2f}%)")
    axs[0].set_ylabel(f"PC2 ({variance[1] * 100:.2f}%)")

    sns.scatterplot(
        data=pc,
        x="PC2",
        y="PC3",
        hue="Sample group",
        style="Replication number",
        ax=axs[1],
    )
    axs[1].set_xlabel(f"PC2 ({variance[1] * 100:.2f}%)")
    axs[1].set_ylabel(f"PC3 ({variance[2] * 100:.2f}%)")
    axs[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.0)

    if plot_ellipses:
        for pc_dims, ax in zip([(1, 2), (2, 3)], axs):
            n = len(pc["Sample group"].unique())
            row_colors = sns.color_palette("tab10", n)
            for group in pc["Sample group"].unique():
                x = pc.loc[pc["Sample group"] == group, f"PC{pc_dims[0]}"]
                y = pc.loc[pc["Sample group"] == group, f"PC{pc_dims[1]}"]
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
                ell.set_facecolor(
                    row_colors[pc["Sample group"].unique().tolist().index(group)]
                )
                ell.set_edgecolor("grey")
                ax.add_artist(ell)

    fig.tight_layout()
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
