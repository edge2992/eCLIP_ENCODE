import pandas as pd
from tqdm import tqdm

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


def _create_accession_count_dict(report: pd.DataFrame, count_method):
    """Accession -> count の辞書を作成する"""
    dd = {}
    for _, row in report.iterrows():
        # Accessionごと
        label = row["Accession"]
        dd[label] = count_method(row)
    return dd


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
