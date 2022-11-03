import pandas as pd
from tqdm import tqdm
from typing import Callable, Dict
from joblib import Parallel, delayed
from src.util.bedfile import read_annotated_bed
from src.util.get_bed_path import get_formatted_file_path
from src.util.bed_format_strategy import FormatStrategy

# python 3.9 is needed
# from collections.abc import Callable


def create_gene_biosample_accession_count_dict(report: pd.DataFrame, count_method):
    """遺伝子 -> 細胞株 -> Accession -> count の辞書を作成する"""
    genes = report["Target label"].unique()
    d = {}
    for gene in tqdm(genes):
        # 遺伝子ごと
        report_ex = report[report["Target label"] == gene]
        biosamples = report_ex["Biosample name"].unique()
        dfs = {}
        for sample in biosamples:
            # 細胞株ごと
            dfs[sample] = _create_accession_count_dict(
                report_ex[report_ex["Biosample name"] == sample], count_method  # type: ignore
            )
        d[gene] = dfs
    return d


def get_geneid_from_assay(row, how):
    df = read_annotated_bed(get_formatted_file_path(row, how))
    return df["gene_id"].dropna().unique().tolist()  # type: ignore


def _create_accession_value(
    report: pd.DataFrame,
    process_method,
    label_method: Callable[[pd.DataFrame], pd.Series] = lambda report: report[
        "Accession"
    ],
):
    """rowを引数とするprocess_methodを適用して、
    アッセイごとになにか値を持つ辞書を作成する
    process_methodは求めたい値を返す関数 (入力: アッセイ情報のrow)
    label_methodはlabelを返す関数 (入力: report)
    """

    # 並列化
    accession_value = Parallel(n_jobs=5, verbose=3)(
        delayed(process_method)(row) for _, row in report.iterrows()
    )

    if accession_value is None:
        raise ValueError("accession_gene is None")

    labels = label_method(report)

    # format to dict
    return {key: value for key, value in zip(labels, accession_value)}


def _create_accession_count_dict(report: pd.DataFrame, count_method):
    """Accession -> count の辞書を作成する"""
    return _create_accession_value(report, count_method)


def convert_gene_biosample_accession_count_dict_to_df(
    gene_biosample_accession_count_dict: dict, report: pd.DataFrame
):
    """遺伝子 -> 細胞株 -> Accession -> count の辞書を
    {key: datasetID -> value : counts in each replicates}のDataFrameに変換する"""
    dd = {}
    report_tags = report[
        ["Dataset", "Target label", "Biosample name"]
    ].drop_duplicates()
    for dataset, gene, biosample in report_tags.itertuples(index=False):
        dd[dataset] = _create_replicate_count_dict(
            gene_biosample_accession_count_dict[gene][biosample],
            report[report["Dataset"] == dataset],
        )
    return pd.DataFrame(dd).T


def _create_replicate_count_dict(accession_count_dict, report):
    """accession -> count の辞書から
    replicates -> count の辞書を作成する"""
    d = {}
    for _, row in report.iterrows():
        label = "replicates_{}".format(row["Biological replicates"])
        d[label] = accession_count_dict[row["Accession"]]
    return d


def count_gene_by_geneType(row: pd.Series, how=FormatStrategy.MAX) -> Dict[str, int]:
    """遺伝子の種類を種類ごとにDict[str, int]で返す"""
    df = read_annotated_bed(get_formatted_file_path(row, how))
    return df["gene_type"].dropna().value_counts().to_dict()  # type: ignore
