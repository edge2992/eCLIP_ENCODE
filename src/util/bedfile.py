import pandas as pd
import os
from dotenv import load_dotenv
from typing import Dict

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

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

COLUMN_INTENDED_ATTRIBUTES = ["gene_id", "gene_name", "gene_type"]


def read_eCLIP_bed(filename: str) -> pd.DataFrame:
    """eCLIPのbedfileにcolumn名を付けて読み込む"""
    return pd.read_table(filename, names=COLUMN_BED_NARROW_PEAK)


def read_intersected_bed(filename: str) -> pd.DataFrame:
    """intersectBedで追加されたcolumnを読み込む"""
    column = COLUMN_BED_NARROW_PEAK + ["gff_" + column for column in COLUMN_ENSEMBL_GFF]
    return pd.read_table(filename, names=column)


def read_annotated_bed(filename: str) -> pd.DataFrame:
    """整形済みの遺伝子情報付きbedfileを読み込む"""
    column = COLUMN_BED_NARROW_PEAK + COLUMN_INTENDED_ATTRIBUTES
    return pd.read_table(filename, names=column)


def count_file_length(filename):
    """ファイルの行数を数える"""
    with open(filename, "r") as f:
        return len(f.readlines())


def count_gene_nunique(filename: str):
    """遺伝子の種類数を数える

    Args:
        df (pd.DataFrame): format_gene_maxで生成されたbedfile

    Returns:
        _type_: _description_
    """
    df = read_annotated_bed(filename).dropna(axis=0, how="any")
    return df["gene_id"].nunique()


def load_report():
    """レポートファイルを読み込んでBiological replicatesでソートする"""
    report = pd.read_table(os.path.join(PROJECT_PATH, "data", "report.tsv"), skiprows=1)  # type: ignore
    report.sort_values("Biological replicates", inplace=True)
    return report


def transform_attribute_to_dict(attribute: str) -> Dict[str, str]:
    """gtfのattributeをdictに変換する"""
    splited_attributes = [x.strip() for x in attribute.replace('"', "").split(";")[:-1]]
    return dict(map(lambda x: x.split(" "), splited_attributes))
