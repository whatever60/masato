import os
import queue
import tempfile
import subprocess
from itertools import cycle

import numpy as np
import pandas as pd
import dendropy
import matplotlib.pyplot as plt
from tqdm.auto import tqdm


def get_coords(
    node: dendropy.Node, max_x: float, overwrite: bool = False
) -> tuple[float, float]:
    if not hasattr(node, "attributes"):
        node.attributes = {}
    attributes = node.attributes
    x = attributes.get("x")
    y = attributes.get("y")
    if not overwrite and x is not None and y is not None:
        return x, y

    x = max_x - node.distance_from_root()

    if node.is_leaf():
        if y is None:
            raise ValueError(
                f"Leaf node {node} has no y value. Set leaf y value first."
            )
    elif node.num_child_nodes() == 0:
        y = -999
        attributes["plot"] = 0
    else:
        y = (
            sum(
                [
                    get_coords(child, max_x=max_x, overwrite=overwrite)[1]
                    for child in node.child_node_iter()
                ]
            )
            / node.num_child_nodes()
        )
    attributes.update({"x": x, "y": y, "plot": 1})
    return x, y


def _get_node_label(node: dendropy.Node) -> str:
    return node.taxon.label if node.taxon else node.label


def get_taxonomy_tree(
    taxon_df: pd.DataFrame, default_branch_length: float = 1.0
) -> dendropy.Tree:
    with (
        tempfile.NamedTemporaryFile(mode="w", suffix=".tsv") as temp_tsv,
        tempfile.TemporaryDirectory() as temp_dir,
    ):
        taxon_df.iloc[:, ::2].apply(lambda x: ";".join(x), axis=1).to_csv(
            temp_tsv.name, header=None, sep="\t"
        )
        proc = subprocess.run(
            [
                "gappa",
                "prepare",
                "taxonomy-tree",
                "--keep-singleton-inner-nodes",
                "--keep-inner-node-names",
                # "--replace-invalid-chars",
                "--allow-file-overwriting",
                "--file-prefix",
                "whatever_",
                "--taxon-list-file",
                temp_tsv.name,
            ],
            cwd=temp_dir,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if proc.returncode != 0:
            raise ValueError("Error in running gappa prepare taxonomy-tree.")

        tree = dendropy.Tree.get(
            path=os.path.join(temp_dir, "whatever_taxonomy_tree.newick"),
            schema="newick",
            rooting="force-rooted",
            preserve_underscores=False,
        ).extract_tree(
            node_filter_fn=lambda n: n.taxon.label in taxon_df.index,
            is_apply_filter_to_leaf_nodes=True,
            suppress_unifurcations=False,
        )
        for edge in tree.postorder_edge_iter():
            if edge.length is None:
                edge.length = default_branch_length

    return tree


def _calc_y(
    tree: dendropy.Tree, taxon_df: pd.DataFrame, collapse_level: str, base: float = 0
) -> dict[str, float]:
    taxon_zotu_count = taxon_df[collapse_level].value_counts().to_dict()
    zotu2taxon = taxon_df[collapse_level].to_dict()
    taxon2y = {}
    taxon2ys = {
        t: np.linspace(-0.4, 0.4, c) if c > 1 else [0]
        for t, c in taxon_zotu_count.items()
    }
    taxon_count_now = {}
    stack = [tree.seed_node]
    y = base
    ordered_taxa = []
    while stack:
        node = stack.pop()
        node_name = _get_node_label(node)
        if node_name in taxon_zotu_count:
            ordered_taxa.append(node_name)
        if node.is_leaf():
            node_taxon = zotu2taxon[node_name]
            if node_taxon not in taxon2y:
                taxon2y[node_taxon] = y
                taxon_count_now[node_taxon] = 0
                y += 1
            node.attributes = {
                "y": taxon2y[node_taxon]
                + taxon2ys[node_taxon][taxon_count_now[node_taxon]]
            }
            taxon_count_now[node_taxon] += 1
        else:
            child_nodes = np.array(node.child_nodes(), dtype=object)
            names = np.array(
                [name for n in child_nodes if (name := _get_node_label(n)) is not None]
            )
            order = np.argsort(names)
            stack.extend(child_nodes[order][::-1].tolist())
    return ordered_taxa


def get_taxon2color(
    taxon_df: pd.DataFrame, levels_of_interest: list[str] = None
) -> dict[str, str]:
    """tab20c for genus, tab20b for family, set3 for order, dark2 for class, set1 for phylum."""
    if levels_of_interest is None:
        levels = ["phylum", "class", "order", "family", "genus"]
    else:
        levels = levels_of_interest
    # colors should be circular in case there are more taxon than colors
    cmaps = ["Set1", "tab20b", "tab20c", "Set3", "Dark2"]
    colors = list(map(lambda c: cycle(plt.get_cmap(c).colors), cmaps))
    level2color = {l: c for l, c in zip(levels, colors)}
    res = {}
    for level, color in level2color.items():
        for taxon in taxon_df[level].unique():
            res[taxon] = next(color)
    res["unknown"] = None
    return res


def get_taxon2marker(
    taxon_df: pd.DataFrame, levels_of_interest: list[str] = None
) -> dict[str, str]:
    # hexagon for domain, pentagon for phylum, square for class, triangle for order, star for family, circle for genus
    markers = "*hps^o"
    if levels_of_interest is None:
        markers = markers
        levels = ["domain", "phylum", "class", "order", "family", "genus"]
    else:
        levels = levels_of_interest
        markers = markers[: len(levels)]
    level2marker = {l: m for l, m in zip(levels, markers)}
    level2marker["species"] = None
    res = {}
    for level, marker in level2marker.items():
        for taxon in taxon_df[level].unique():
            res[taxon] = marker
    res["unknown"] = None
    return res


def plot_tree(
    tree: dendropy.Tree,
    node_label2marker: dict[str, float],
    node_label2color: dict[str, str],
    node_label2alpha: dict[str, float],
    ax: plt.Axes,
    color_propagate: bool = False,
    annotate: bool = False,
    nodes_to_drop: list[str] = tuple(),
) -> None:
    """Plot dendrogram with matplotlib.
    Plotting starts from root to leaf, and each iteration draws a horizontal line from parent
        coordinates to the x coordinate of its children, and a vertical line that spans
        the y coordinates of its children. Along the way, make sure that
        each child has the same x coordinate, otherwise the tree is ill-formed.
    """
    node_label2color = dict(node_label2color)
    # nodes_to_plot = queue.Queue()
    # nodes_to_plot.put(tree.seed_node)
    # while not nodes_to_plot.empty():
        # node = nodes_to_plot.get()
    nodes_to_plot = [tree.seed_node]
    while nodes_to_plot:
        node = nodes_to_plot.pop()
        node_name = _get_node_label(node)
        if not node.attributes.get("plot"):
            continue
        x, y = get_coords(node, max_x=None)

        # if node.parent_node is not None:
        #     xp, yp = get_coords(node.parent_node, max_x=None)
        #     # horizontal line from parent to child
        #     ax.plot(
        #         [xp, x],
        #         [y, y],
        #         "k",
        #         alpha=edge_label2alpha.get(node.edge.label, 1),
        #         solid_capstyle="round",
        #         lw=1,
        #     )
        #     if x == 0:
        #         plotted_nodes.append(node)
        # x coord of child nodes

        color = node_label2color.get(node_name, "k")
        if node.num_child_nodes():
            child_coords = np.array(
                [
                    get_coords(child, max_x=None)
                    for child in node.child_node_iter()
                    if child.attributes.get("plot")
                ]
            )
            order = child_coords[:, 1].argsort()
            child_nodes = np.array(node.child_nodes(), dtype=object)[order]
            child_coords = child_coords[order]

            if len(child_coords) > 1:
                # plot a vertical line from the y coord of the first child to the last child
                if len(np.unique(child_coords[:, 0])) > 1:
                    # validate that all child nodes have the same x coord
                    raise ValueError(
                        f"Child nodes of {node} have different x coordinates."
                    )
            for i, (node_child, (xc, yc)) in enumerate(zip(child_nodes, child_coords)):
                child_name = _get_node_label(node_child)
                if not node_child.attributes.get("plot") or child_name in nodes_to_drop:
                    continue
                yt = (
                    max(y, child_coords[i - 1, 1])
                    if yc > y or len(child_coords) == 1
                    else min(y, child_coords[i + 1, 1])
                )

                # alpha goes with child node, color goes with parent node
                ax.plot(
                    [x, x, xc],
                    [yt, yc, yc],
                    alpha=node_label2alpha.get(child_name, 1),
                    color=color,
                    lw=1,
                    solid_capstyle="round",
                    solid_joinstyle="round",
                )
                # if node_name == "unknown":
                #     import pdb; pdb.set_trace()
                if color_propagate and (
                    not node_label2color.get(child_name)
                    or not node_label2marker.get(child_name)
                ):
                    node_label2color[child_name] = color
                nodes_to_plot.append(node_child)
        marker = node_label2marker.get(node_name)
        if marker is not None:
            ax.scatter(
                x,
                y,
                s=5,
                color=color,
                alpha=node_label2alpha.get(node_name, 1),
                marker=marker,
                zorder=10,
            )
        if annotate:
            ax.text(
                x,
                y,
                node_name,
                fontsize=3,
                ha="right",
                va="center",
                color="k",
                zorder=10,
            )
