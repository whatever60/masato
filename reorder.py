import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs


def hierarchical_order(
    df: pd.DataFrame, num_clusters: int, method: str = "ward"
) -> (np.ndarray, np.ndarray):
    """
    Perform hierarchical clustering on the provided DataFrame and return the ordered indices and cluster labels.

    Parameters:
    df (pd.DataFrame): The DataFrame on which to perform hierarchical clustering.
    num_clusters (int): The number of clusters to divide the data into.
    method (str): The linkage method to use for clustering.

    Returns:
    tuple: A tuple containing:
           - np.ndarray: Array of ordered indices based on hierarchical clustering.
           - np.ndarray: Array of cluster labels for each row in the DataFrame.
    """
    linkage_result = linkage(df, method=method)
    ordered_indices = leaves_list(linkage_result)
    cluster_labels = (
        fcluster(linkage_result, num_clusters, criterion="maxclust")
        if num_clusters is not None
        else None
    )
    return ordered_indices, cluster_labels


def get_order(
    df: pd.DataFrame,
    num_clusters: int,
    inter_hier: bool = None,
    intra_hier: bool = None,
    seed: int = 42,
    method: str = "kmeans",
) -> (np.ndarray, np.ndarray):
    """
    Perform clustering and reorder rows of the dataframe based on hierarchical clustering or KMeans as specified.

    Parameters:
    df (pd.DataFrame): The dataframe to be clustered and reordered.
    num_clusters (int): The number of clusters to form.
    inter_hier (bool): Whether to reorder clusters based on hierarchical clustering of centroids for KMeans.
    intra_hier (bool): Whether to reorder within clusters based on hierarchical clustering for KMeans.
    seed (int): Random seed for reproducibility.
    method (str): Clustering method, 'kmeans' or 'hier' for hierarchical.

    Returns:
    tuple: A tuple containing:
           - np.ndarray: Array of row indices that reorder the dataframe rows.
           - np.ndarray: Array of cluster labels for each row in the input dataframe.
    """
    df = df.reset_index(drop=True)
    df_filled = df.fillna(0)

    if method == "kmeans":
        kmeans = KMeans(n_clusters=num_clusters, random_state=seed)
        cluster_labels = kmeans.fit_predict(df_filled)

        if inter_hier:
            centroids = kmeans.cluster_centers_
            ordered_indices, _ = hierarchical_order(pd.DataFrame(centroids), None)
            label_map = {
                old_label: new_label
                for new_label, old_label in enumerate(ordered_indices)
            }
            ordered_labels = np.array([label_map[label] for label in cluster_labels])
        else:
            ordered_labels = cluster_labels

        ordered_df = pd.DataFrame(df_filled).assign(label=ordered_labels)
        ordered_df.sort_values("label", inplace=True)
        new_indices = []
        if intra_hier:
            for label in np.unique(ordered_labels):
                cluster_indices = ordered_df[ordered_df["label"] == label].index
                cluster_data = df_filled.loc[cluster_indices]
                if len(cluster_data) == 1:
                    new_indices.extend(cluster_indices)
                else:
                    intra_cluster_order, _ = hierarchical_order(cluster_data, None)
                    new_indices.extend(cluster_indices[intra_cluster_order].tolist())
            new_indices = np.array(new_indices)
        else:
            new_indices = ordered_df.index.to_numpy()

    elif method == "hier":
        ordered_indices, cluster_labels = hierarchical_order(df_filled, num_clusters)
        ordered_df = df_filled.iloc[ordered_indices]
        new_indices = ordered_df.index.to_numpy()

    return np.array(new_indices), cluster_labels


if __name__ == "__main__":
    # Testing with simulated 2D Gaussian blobs
    X, y = make_blobs(n_samples=200, centers=4, cluster_std=0.60, random_state=0)
    blob_df = pd.DataFrame(X, columns=["Feature1", "Feature2"])

    # Apply the function
    indices, labels = get_order(
        blob_df, num_clusters=4, inter_hier=True, intra_hier=True, method="kmeans"
    )
    print(indices.shape, labels.shape)
