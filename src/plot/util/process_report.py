# レポートを処理して、描画しやすいデータを返す関数を定義する
from typing import Callable
import pandas as pd
from src.util.bed_format_strategy import FormatStrategy
from src.plot.util.process_by_accession import (
    _create_accession_value,
    count_gene_by_geneType,
    get_geneid_from_assay,
)


def label_full(report: pd.DataFrame):
    return (
        report["Target label"]
        + " "
        + report["Biosample name"]
        + " "
        + report["Accession"]
    )


def label_protein_biosample(report: pd.DataFrame) -> pd.Series:
    return report["Target label"] + " " + report["Biosample name"]


def count_gene(
    report: pd.DataFrame, label_method: Callable[[pd.DataFrame], pd.Series] = label_full
) -> pd.DataFrame:
    """遺伝子の種類を種類ごとに数える"""
    """行がACCESSIONで, 列がgene_typeのDataFrameを返す"""
    # 並列処理
    result = _create_accession_value(
        report,
        lambda row: count_gene_by_geneType(row, FormatStrategy.MAX),
        label_method,
    )

    df = pd.DataFrame(result)
    return df.T


def get_gene_ids(
    report: pd.DataFrame,
    label_method: Callable[[pd.DataFrame], pd.Series] = lambda x: x["Accession"],
):
    """
    アッセイごとに遺伝子のidのListを取得する
    アッセイごとにユニークなgene_idのリストを調べる
    Dict[accession, List[gene]を作成する
    """
    return _create_accession_value(
        report, lambda row: get_geneid_from_assay(row, FormatStrategy.MAX), label_method
    )
