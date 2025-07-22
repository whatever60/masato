import warnings

import anndata as ad
import biom
import pandas as pd

from scipy import sparse


def should_warn_metadata_loss(obj) -> bool:
    """Determine if losing metadata would actually cause data loss.

    Args:
        obj: AnnData, biom.Table, or similar object.

    Returns:
        True if metadata is non-empty and a warning should be issued.
    """
    if isinstance(obj, ad.AnnData):
        obs_has_data = not obj.obs.empty
        var_has_data = not obj.var.empty
        return obs_has_data or var_has_data

    if isinstance(obj, biom.Table):
        obs_has_metadata = sample_has_metadata = False
        meta_obs = obj.metadata(axis="observation")
        meta_sample = obj.metadata(axis="sample")
        if isinstance(meta_obs, tuple):
            obs_has_metadata = any(table_md is not None for table_md in meta_obs)
        if isinstance(meta_sample, tuple):
            sample_has_metadata = any(table_md is not None for table_md in meta_sample)
        return obs_has_metadata or sample_has_metadata

    return False


def read_table(
    table_path: str,
    index_col: str | int = 0,  # Only used for reading flat files
    comment: str | None = None,  # Only used for reading flat files
    dtype: str = "int",  # If return_type is biom, dtype is forced to be float.
    index_name: str | None = None,
    return_type: str = "df",
) -> pd.DataFrame | ad.AnnData | biom.Table:
    """Read a table from a file and return it as a DataFrame, AnnData, or biom.Table.

    Args:
        table_path: Path to the table file. Could be a tsv, csv or biom.
        index_col: Column to use as the row labels of the DataFrame.
        comment: Character to split comments in flat files.
        dtype: Data type to force.
        index_name: Name to set for index.
        return_type: "df" for DataFrame, "adata" for AnnData, "biom" for biom.Table.

    Returns:
        DataFrame, AnnData or biom.Table depending on return_type.
    """
    if table_path.endswith(".biom"):
        table_biom = biom.load_table(table_path)

        if return_type == "biom":  # biom -> biom (sparse)
            return table_biom

        if return_type == "adata":  # biom -> AnnData (sparse)
            # if should_warn_metadata_loss(table_biom):
            #     warnings.warn(
            #         f"Metadata lost when reading {table_path} into AnnData: only IDs preserved."
            #     )
            try:
                obs = table_biom.metadata_to_dataframe(axis="sample")
            except KeyError:
                obs = pd.DataFrame(index=table_biom.ids(axis="sample"))
            try:
                var = table_biom.metadata_to_dataframe(axis="observation")
            except KeyError:
                var = pd.DataFrame(index=table_biom.ids(axis="observation"))

            matrix = table_biom.matrix_data.T.astype(dtype)
            adata = ad.AnnData(X=matrix, obs=obs, var=var)
            adata.var.index.name = (
                table_biom.table_id if index_name is None else index_name
            )
            return adata

        if return_type == "df":  # biom -> DataFrame (dense)
            if should_warn_metadata_loss(table_biom):
                warnings.warn(
                    f"Metadata lost when reading {table_path} into DataFrame: only IDs preserved."
                )
            df = table_biom.to_dataframe(dense=True).astype(dtype)
            df.index.name = table_biom.table_id if index_name is None else index_name
            return df

        raise ValueError("return_type must be one of: 'df', 'adata', 'biom'.")

    else:  # flat file (return data is always dense)
        if table_path.endswith(".csv") or table_path.endswith(".csv.gz"):
            sep = ","
        elif table_path.endswith(".tsv") or table_path.endswith(".tsv.gz"):
            sep = "\t"
        else:
            raise ValueError("Unsupported table file format.")

        df = pd.read_csv(
            table_path,
            sep=sep,
            index_col=index_col,
            comment=comment,
            dtype={index_col: str},
        )
        if index_name is not None:
            df.index.name = index_name

        if return_type == "df":  # flat file -> DataFrame
            return df
        elif return_type == "adata":  # flat file -> AnnData
            # warnings.warn(
            #     f"Flat file {table_path} has no metadata, creating AnnData from matrix only."
            # )
            return ad.AnnData(
                X=df.values.T.astype(dtype),
                obs=pd.DataFrame(index=df.columns),
                var=pd.DataFrame(index=df.index),
            )
        elif return_type == "biom":  # flat file -> biom
            # warnings.warn(
            #     f"Flat file {table_path} has no metadata, creating biom.Table from matrix only."
            # )
            return biom.Table(
                df.values.astype(dtype),
                observation_ids=df.index.tolist(),
                sample_ids=df.columns.tolist(),
                table_id=df.index.name,
                type=dtype,
            )
        else:
            raise ValueError("return_type must be one of: 'df', 'adata', 'biom'.")


def write_table(table: pd.DataFrame | ad.AnnData | biom.Table, table_path: str) -> None:
    """Write a table to a file.

    Args:
        table: DataFrame, AnnData, or biom.Table to be written.
        table_path: Path to the output file. Could be a tsv, csv or biom.
    """
    if table_path.endswith(".biom"):
        if isinstance(table, biom.Table):  # always float
            data = table
        elif isinstance(table, ad.AnnData):
            # if should_warn_metadata_loss(table):
            #     warnings.warn(
            #         "Writing AnnData to biom format will lose obs/var metadata."
            # )
            matrix = table.X.T
            data = biom.Table(
                matrix,
                observation_ids=table.var_names.to_list(),
                sample_ids=table.obs_names.to_list(),
                observation_metadata=table.var.values.tolist(),
                sample_metadata=table.obs.values.tolist(),
                table_id=table.var.index.name,
            )
        elif isinstance(table, pd.DataFrame):
            data = biom.Table(
                table.values,
                observation_ids=table.index.tolist(),
                sample_ids=table.columns.tolist(),
                table_id=table.index.name,
            )
        else:
            raise ValueError("Unsupported input type for biom output.")
        with biom.util.biom_open(table_path, "w") as f:
            data.to_hdf5(f, "whatever60")

    else:  # flat file
        if table_path.endswith(".csv") or table_path.endswith(".csv.gz"):
            sep = ","
        elif table_path.endswith(".tsv") or table_path.endswith(".tsv.gz"):
            sep = "\t"
        else:
            raise ValueError("Unsupported table file format.")

        if isinstance(table, biom.Table):  # always float
            if should_warn_metadata_loss(table):
                warnings.warn(
                    f"Writing biom.Table to flat file {table_path} will lose metadata, saving only matrix."
                )
            df = table.to_dataframe(dense=True)
        elif isinstance(table, ad.AnnData):
            if should_warn_metadata_loss(table):
                warnings.warn(
                    f"Writing AnnData to flat file {table_path} will lose obs/var metadata, saving only X."
                )
            df = pd.DataFrame(
                table.X.toarray().T if sparse.issparse(table.X) else table.X.T,
                index=table.var_names,
                columns=table.obs_names,
            )
        elif isinstance(table, pd.DataFrame):
            df = table
        else:
            raise ValueError("Unsupported input type for csv/tsv output.")

        df.to_csv(table_path, sep=sep)
