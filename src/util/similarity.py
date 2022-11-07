import os
from typing import Callable, Dict, List
import pandas as pd
from dotenv import load_dotenv
from src.plot.util.process_intersect_gene import count_interection
from src.plot.util.process_report import (
    count_gene,
    gene_ids_eCLIP,
    label_protein_biosample,
)
from src.util.uniprot import idmapping

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
CLUSTALO_DISTMAT_PATH = os.environ["CLUSTALO_DISTMAT_PATH"]


def similarity_set(gene_list1, gene_list2):
    """二つの遺伝子の集合の類似度を計算する"""
    # TODO: 実装方針を立てる
    return len(set(gene_list1) & set(gene_list2))


def load_protein_sequence_similarity():
    # clustal Omegaのpercent identity Matrixを読み取って、pd.DataFrameにする
    similarities: List[List[str]] = []
    labels: List[str] = []
    mapping_dict: Dict[str, str] = {value: key for key, value in idmapping().items()}
    with open(CLUSTALO_DISTMAT_PATH, "r") as f:
        lines = f.readlines()
        for line in lines[1:]:
            if line.startswith("#"):
                continue
            elif line == "\n":
                continue
            else:
                similarities.append(line.split()[1:])
                labels.append(line.split()[0])
    df = (
        pd.DataFrame(similarities, index=labels, columns=labels, dtype=float)
        .loc[list(mapping_dict.keys()), list(mapping_dict.keys())]
        .rename(index=mapping_dict, columns=mapping_dict)
    )
    return df


def lift_protein(
    report: pd.DataFrame,
    label_method: Callable[[pd.DataFrame], pd.Series] = label_protein_biosample,
) -> pd.DataFrame:
    """タンパク質のリフト値を計算する。
    リフト値はマーケット分析で使われる指標
    eCLIPはタンパク質を中心に解析をするが、
    視点を変えてRNAごとに結合するタンパク質が出現する頻度を考えて、
    タンパク質同士の出現の相関を求める。
    出現の仕方が互いに一切関係なく、独立な場合にはリフト値は1となる。
    正の相関だと1より大きくなり負の相関だと1より小さくなる。
    """
    intersect = count_interection(report, label_method)
    gene = count_gene(report, label_method)
    all_gene_count = len(gene_ids_eCLIP(report))
    columns = intersect.columns

    gene_value = gene[columns].to_numpy()
    value = (
        intersect.to_numpy()
        / (gene_value * gene_value.reshape((-1, 1)))
        * all_gene_count
    )
    return pd.DataFrame(value, index=columns, columns=columns)
