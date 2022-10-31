import pandas as pd
import os


PROJECT_PATH = "/mnt/H/MYWORK/eCLIP_ENCODE"

COLUMN_BED_NARROW_PEAK = [
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

COLUMN_ENSEMBL_GFF = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]


def read_eCLIP_bed(filename: str) -> pd.DataFrame:
    """eCLIPのbedfileにcolumn名を付けて読み込む"""
    return pd.read_table(filename, names=COLUMN_BED_NARROW_PEAK)


def read_intersected_bed(filename: str) -> pd.DataFrame:
    """intersectBedで追加されたcolumnを読み込む"""
    column = COLUMN_BED_NARROW_PEAK + ["gff_" + column for column in COLUMN_ENSEMBL_GFF]
    return pd.read_table(filename, names=column)


def count_file_length(filename):
    """ファイルの行数を数える"""
    with open(filename, "r") as f:
        return len(f.readlines())


def load_report():
    """レポートファイルを読み込んでBiological replicatesでソートする"""
    report = pd.read_table(os.path.join(PROJECT_PATH, "data", "report.tsv"), skiprows=1)
    report.sort_values("Biological replicates", inplace=True)
    return report
