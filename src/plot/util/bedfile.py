import pandas as pd
import os

PROJECT_PATH = "/mnt/H/MYWORK/eCLIP_ENCODE"


def read_eCLIP_bed(filename):
    """eCLIPのbedfileにcolumn名を付けて読み込む"""
    columns = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "singleValue",
        "pValue",
        "qValue",
        "peak",
    ]
    df = pd.read_table(filename, names=columns)
    return df


def count_file_length(filename):
    """ファイルの行数を数える"""
    with open(filename, "r") as f:
        return len(f.readlines())


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


def load_report():
    """レポートファイルを読み込んでBiological replicatesでソートする"""
    report = pd.read_table(os.path.join(PROJECT_PATH, "data", "report.tsv"), skiprows=1)
    report.sort_values("Biological replicates", inplace=True)
    return report
