from enum import Enum
import pandas as pd


# use this for format_gene_binding_sites
FormatStrategy = Enum("FormatStrategy", ["MAX"])


def format_max(intersected: pd.DataFrame):
    """singleValueがmaxのbindingsiteを抽出して整形する"""
    assert intersected["singleValue"].dtype == "float64"
    return (
        intersected.sort_values("singleValue", ascending=False)
        .drop_duplicates(["gff_seqname", "gff_start", "gff_end"], keep="first")
        .reset_index()
    )
