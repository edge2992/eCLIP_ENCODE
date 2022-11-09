from src.util.bedfile import read_annotated_bed
from src.util.get_bed_path import get_formatted_file_path
from src.util.bed_format_strategy import FormatStrategy
import pandas as pd


def get_gene_list(report: pd.DataFrame, accession: str):
    """Accessionから遺伝子のリストを取得する"""
    seri = report[report["Accession"] == accession]
    assert len(seri) == 1
    return (
        read_annotated_bed(get_formatted_file_path(seri.iloc[0], FormatStrategy.MAX))[
            "gene_id"
        ]
        .dropna()
        .unique()
        .tolist()
    )
