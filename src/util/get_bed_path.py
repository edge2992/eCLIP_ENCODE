import os
import pandas as pd
from dotenv import load_dotenv

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]


def get_file_path(row: pd.Series):
    """report.tsvの行からファイルのパスを取得する"""
    return os.path.join(
        PROJECT_PATH,
        "data",
        row["Assay term name"],
        row["Target label"],
        row["Biosample name"].split()[0],  # adrenal gland, K562, HepG2
        row["Accession"] + ".bed",
    )


def get_annotated_file_path(row: pd.Series):
    """report.tsvの行からgffとintersectしたファイルのパスを取得する"""
    return os.path.join(
        PROJECT_PATH,
        "annotated_data",
        "intersect",
        row["Assay term name"],
        row["Target label"],
        row["Biosample name"].split()[0],  # adrenal gland, K562, HepG2
        row["Accession"] + ".bed",
    )


def get_formatted_file_path(row: pd.Series, how):
    """report.tsvの行から遺伝子ごとに整形済みのファイルのパスを取得する"""
    return os.path.join(
        PROJECT_PATH,
        "annotated_data",
        "gene",
        how.name.lower(),
        row["Assay term name"],
        row["Target label"],
        row["Biosample name"].split()[0],  # adrenal gland, K562, HepG2
        row["Accession"] + ".bed",
    )
