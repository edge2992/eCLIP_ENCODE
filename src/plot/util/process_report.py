# レポートを処理して、描画しやすいデータを返す関数を定義する
import pandas as pd
from src.util.bed_format_strategy import FormatStrategy
from src.plot.util.process_by_accession import (
    _create_accession_value,
    count_gene_by_geneType,
)


def count_gene(report: pd.DataFrame):
    """遺伝子の種類を種類ごとに数える"""
    """行がACCESSIONで, 列がgene_typeのDataFrameを返す"""
    # 並列処理
    result = _create_accession_value(
        report,
        lambda row: count_gene_by_geneType(row, FormatStrategy.MAX),
        lambda x: (
            x["Target label"] + " " + x["Biosample name"] + " " + x["Accession"]
        ),  # type: ignore
    )

    df = pd.DataFrame(result)
    return df.T
