# このフォルダで使用する便利関数
import pandas as pd

from src.eclip import Compare, Dataset
from src.plot.util.process_report import count_gene
from src.util.bed_format_strategy import FormatStrategy
from src.util.bedfile import load_replicateIDR_report, read_annotated_bed
from src.util.get_bed_path import get_formatted_file_path
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

    def intersection_keyword(row: pd.Series):
        keyword1 = get_keyword(row["Dataset_1"])
        keyword2 = get_keyword(row["Dataset_2"])
        return list(set(keyword1) & set(keyword2))

    keywords = data.apply(intersection_keyword, axis=1)  # type: ignore
    keywords.name = "keyword"
    return dict(
        keywords.explode().reset_index().groupby("keyword")["index"].apply(list)
    )


def describe_dataset_pair(row: pd.Series):
    dataset1 = Dataset(row["Dataset_1"])
    dataset2 = Dataset(row["Dataset_2"])
    compare = Compare([dataset1, dataset2])

    print("-" * 20)
    print(compare)
    print("gene1: {}".format(len(dataset1.genes)))
    print("gene2: {}".format(len(dataset2.genes)))
    print("gene intersection: {}".format(len(compare.gene_intersection)))
    print("keyword1: {}".format(len(dataset1.keywords)))
    print("keyword2: {}".format(len(dataset2.keywords)))
    print("keyword intersection: {}".format(compare.keyword_intersection))
    print("simpson: {}".format(row["Simpson"]))
