from typing import Callable, List, Dict
import pandas as pd
from tqdm import tqdm

from src.plot.util.process_report import (
    count_gene_type,
    get_gene_ids,
    label_protein_biosample,
)


def count_interection(
    report: pd.DataFrame,
    label_method: Callable[[pd.DataFrame], pd.Series] = label_protein_biosample,
) -> pd.DataFrame:
    """_summary_reportを受け取り、アッセイ間の遺伝子のかぶりを数える

    Args:
        report (pd.DataFrame): _description_
        label_method (Callable[[pd.DataFrame], pd.Series], optional): _description_. Defaults to label_protein_biosample.

    Returns:
        pd.DataFrame: _description_
    """

    accession_genes: Dict[str, List[str]] = get_gene_ids(report, label_method)
    columns = sorted(list(accession_genes.keys()))
    data = pd.DataFrame(columns=columns, index=columns, dtype=int)
    for (k1, v1) in tqdm(accession_genes.items()):
        for (k2, v2) in accession_genes.items():
            data.loc[k1, k2] = len(set(v1) & set(v2))
    return data


def assay_with_many_gene(
    report: pd.DataFrame,
    threshold: int,
    prepare_label: Callable[[pd.DataFrame], pd.Series] = lambda x: x["Accession"],
) -> List[str]:
    """遺伝子数が多いアッセイを取得する"""
    data: pd.Series[int] = count_gene_type(report, prepare_label).sum(axis=1)
    intended_indexes: List[str] = sorted(list(data[data > threshold].index))
    print("{} -> {}".format(data.shape[0], len(intended_indexes)))
    return intended_indexes
