# このフォルダで使用する便利関数
from typing import Dict, List
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from tqdm import tqdm

from src.plot.util.process_report import count_gene
from src.util.bedfile import load_replicateIDR_report, read_annotated_bed
from src.util.similarity_protein import InteractionSimilarity, ProteinSimilarity
from src.util.similarity_strategy import (
    TAPE,
    BlastP,
    Cosine,
    KeywordCosine,
    Lift,
    Simpson,
)
from src.util.get_bed_path import get_formatted_file_path
from src.util.bed_format_strategy import FormatStrategy
from src.util.uniprot import load_keyword_report


def target_report(threshold_gene_num: int, biosample: str):
    """調査対象となる実験のdatasetを用意する
    biosample: 細胞株 HepG2 or K562
    threshold_gene_num: 相互作用する遺伝子の数がこの値以上
    """
    report = load_replicateIDR_report()
    report = report[report["Biosample name"] == biosample].reset_index(drop=True)
    target_report = report[
        (count_gene(report, lambda row: row["Dataset"]) >= threshold_gene_num).to_list()
    ]
    return target_report


def metrics(report: pd.DataFrame):
    """統計値を整形する
    タンパク質の類似度: TAPE(cosine distance), blastp(bit score), keyword(cosine distance)
    相互作用の遺伝子の類似度: simpson index
    """
    p_similairty = ProteinSimilarity()
    inter_similarity = InteractionSimilarity()

    def protein_similarity(strategy):
        p_similairty.setStrategy(strategy)
        return p_similairty.flatten_tri(
            p_similairty.executeStrategy(transform=True), False
        )

    def interaction_similarity(strategy):
        inter_similarity.setStrategy(strategy)
        return inter_similarity.flatten_tri(inter_similarity.executeStrategy(), False)

    REPORT_COLUMNS = ["Dataset", "Target label", "Biosample name"]

    data = pd.DataFrame(
        {
            "TAPE": protein_similarity(TAPE(report=report)),
            "keyword": protein_similarity(KeywordCosine(report=report)),
            "blastp": protein_similarity(BlastP(report=report)),
            "simpson": interaction_similarity(Simpson(report=report)),
            "lift": interaction_similarity(Lift(report=report)),
            "cosine": interaction_similarity(Cosine(report=report)),
        }
    )
    index_n = np.where(np.triu(squareform(np.ones(data.shape[0]))))
    desc = pd.concat(
        [
            report.iloc[index_n[0]]
            .loc[:, REPORT_COLUMNS]
            .reset_index(drop=True)
            .add_suffix("_1"),
            report.iloc[index_n[1]]
            .loc[:, REPORT_COLUMNS]
            .reset_index(drop=True)
            .add_suffix("_2"),
        ],
        axis=1,
    )
    return pd.concat([desc, data], axis=1).sort_values("TAPE").reset_index(drop=True)


def get_geneset(dataset: str, how=FormatStrategy.MAX):
    """タンパク質に結合する遺伝子のセットを取得する"""
    report = load_replicateIDR_report().set_index("Dataset")
    df = read_annotated_bed(get_formatted_file_path(report.loc[dataset], how))
    return list(set(df["gene_name"]))


def get_keyword(dataset: str):
    """タンパク質のキーワードを取得する"""
    keywords = load_keyword_report().set_index("From")
    report = load_replicateIDR_report().set_index("Dataset")
    target: str = report.loc[dataset]["Target label"]
    data_str: str = keywords.loc[target, "Keywords"]  # type: ignore
    return [key.strip() for key in data_str.split(";")]


def convert_to_dict_exp_pair_by_keyword(data: pd.DataFrame):
    """keyword -> dataのindexの辞書を作成する"""
    keyword_experiment_pair: Dict[str, List] = {}

    for index, row in tqdm(data.iterrows()):
        keyword1 = get_keyword(row["Dataset_1"])
        keyword2 = get_keyword(row["Dataset_2"])
        intersection_keyword = list(set(keyword1) & set(keyword2))
        for k in intersection_keyword:
            if k not in keyword_experiment_pair:
                keyword_experiment_pair[k] = []
            keyword_experiment_pair[k].append(index)
    return keyword_experiment_pair


def describe_dataset_pair(row: pd.Series):
    gene1 = get_geneset(row["Dataset_1"])
    gene2 = get_geneset(row["Dataset_2"])
    keyword1 = get_keyword(row["Dataset_1"])
    keyword2 = get_keyword(row["Dataset_2"])

    print("-" * 20)
    print("dataset1: {}".format(row["Dataset_1"]))
    print("dataset2: {}".format(row["Dataset_2"]))
    print("protein1: {}".format(row["Target label_1"]))
    print("protein2: {}".format(row["Target label_2"]))
    print("gene1: {}".format(len(gene1)))
    print("gene2: {}".format(len(gene2)))
    print("keyword1: {}".format(len(keyword1)))
    print("keyword2: {}".format(len(keyword2)))
    print("keyword intersection: {}".format(list(set(keyword1) & set(keyword2))))
    print("intersection: {}".format(len(set(gene1) & set(gene2))))
    print("simpson: {}".format(row["simpson"]))
