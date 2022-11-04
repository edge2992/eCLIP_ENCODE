# レポートを処理して、描画しやすいデータを返す関数を定義する
from typing import Callable
import pandas as pd
from src.util.bed_format_strategy import FormatStrategy
from src.plot.util.process_by_accession import (
    _create_accession_value,
    count_gene_by_geneType,
    get_geneid_from_assay,
)
from src.util.bedfile import read_annotated_bed
from src.util.get_bed_path import get_formatted_file_path


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


def count_gene_type(
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


def count_gene(
    report: pd.DataFrame, label_method: Callable[[pd.DataFrame], pd.Series] = label_full
) -> pd.Series:
    """遺伝子の種類を数える"""
    """行がACCESSIONで, 列がgene_typeのDataFrameを返す"""

    def _count_gene(row: pd.Series, how=FormatStrategy.MAX) -> int:
        """遺伝子の数を返す"""
        df = read_annotated_bed(get_formatted_file_path(row, how))
        return len(df["gene_id"].dropna().unique().tolist())  # type: ignore

    # 並列処理
    result = _create_accession_value(
        report,
        lambda row: _count_gene(row, FormatStrategy.MAX),
        label_method,
    )
    return pd.Series(result)


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
