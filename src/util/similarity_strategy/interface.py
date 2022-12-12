import os
from abc import ABC, abstractmethod
from itertools import chain
from typing import Callable, Dict, List, Union

import numpy as np
import pandas as pd
from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]


class SimilarityStrategy(ABC):
    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
        label_method: Callable[[pd.DataFrame], pd.Series] = lambda df: df["Dataset"],
    ):
        if report is None:
            from src.util.bedfile import load_replicateIDR_report

            self.report = load_replicateIDR_report()
        else:
            self.report = report
        self.loadfile = loadfile
        self.label_method = label_method

    @abstractmethod
    def execute(self) -> pd.DataFrame:
        raise NotImplementedError()


class ProteinSimilarityStrategy(SimilarityStrategy):

    UNIPROT_TABLE = os.path.join(PROJECT_PATH, "data/uniprot", "reviewed.tsv")

    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
        label_method: Callable[[pd.DataFrame], pd.Series] = lambda df: df["Dataset"],
    ):
        super().__init__(report, loadfile, label_method)

    def execute(self) -> pd.DataFrame:
        similarities = self._protein_similarity()
        mapping_dict = self._idmapping()

        return (
            similarities.loc[list(mapping_dict.keys()), list(mapping_dict.keys())]
            .rename(index=mapping_dict, columns=mapping_dict)
            .sort_index(axis=0)
            .sort_index(axis=1)
        )

    def transform(self, similarities: pd.DataFrame) -> pd.DataFrame:
        """reportに沿って、タンパク質の類似度ベクトルを作成する"""
        from scipy.spatial.distance import squareform

        report_dataset_indexed = self.report.set_index("Dataset")

        def get_protein_from_dataset(dataset: str) -> str:
            """データセットからタンパク質名を取得する"""
            seri = report_dataset_indexed.loc[dataset]
            return str(seri["Target label"])

        def get_similarity(dataset1: str, dataset2: str):
            gene1 = get_protein_from_dataset(dataset1)
            gene2 = get_protein_from_dataset(dataset2)
            return similarities.loc[gene1, gene2]

        def upper_triangle(label: List[str]):
            # Accessionのラベルにしたがって、
            # タンパク質の類似度行列を並び替える
            # scipy _pdist_callable like function
            # https://github.com/scipy/scipy/blob/v1.9.3/scipy/spatial/distance.py#L1943-L2246
            n = len(label)
            out_size = (n * (n - 1)) // 2
            dm = np.empty(out_size, dtype=np.float32)
            k = 0
            for i in range(n - 1):
                for j in range(i + 1, n):
                    dm[k] = get_similarity(label[i], label[j])
                    k += 1
            return dm

        labels: List[str] = self.report["Dataset"].to_list()
        return pd.DataFrame(
            squareform(upper_triangle(labels)), index=labels, columns=labels
        )

    @abstractmethod
    def _protein_similarity(self) -> pd.DataFrame:
        # DataFrameとColumnとIndexの順番は自由
        raise NotImplementedError()

    def _idmapping(self) -> Dict[str, str]:
        """uniprot idmappingからダウンロードしたテーブルデータを使用して、
        uniprotのproteinidとeCLIPで使われているタンパク質の対応表を作成する
        example: sp|O00425|IF2B3_HUMAN -> IGF2BP3
        """

        df = (
            pd.read_table(self.UNIPROT_TABLE)
            .sort_values("Length", ascending=False)
            .drop_duplicates("From")
        )
        df["label"] = "sp|" + df["Entry"] + "|" + df["Entry Name"]
        return df[["From", "label"]].set_index("label").to_dict()["From"]  # type: ignore


class InteractionSimilarityStrategy(SimilarityStrategy):
    _interaction_intersection: Union[None, pd.DataFrame] = None
    _interaction_union: Union[None, pd.DataFrame] = None
    _accession_genes: Union[None, Dict[str, List[str]]] = None
    _nunique_gene: Union[None, int] = None

    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
        label_method: Callable[[pd.DataFrame], pd.Series] = lambda df: df["Dataset"],
    ):
        super().__init__(report, loadfile, label_method)

    @property
    def accession_genes(self) -> Dict[str, List[str]]:
        # 相互作用数と遺伝子数をそれぞれ読み込む

        def cached_load(
            report: pd.DataFrame,
        ) -> Dict[str, List[str]]:
            from src.plot.util.process_report import get_gene_ids

            accession_genes: Dict[str, List[str]] = get_gene_ids(
                report, self.label_method
            )
            return accession_genes

        if not isinstance(self._accession_genes, dict):
            self._accession_genes = cached_load(self.report)

        return self._accession_genes

    @property
    def interaction_intersection(self) -> pd.DataFrame:
        if self._interaction_intersection is None:
            data = pd.DataFrame(
                columns=list(self.accession_genes.keys()),
                index=list(self.accession_genes.keys()),
                dtype=int,
            )
            for (k1, v1) in tqdm(self.accession_genes.items()):
                for (k2, v2) in self.accession_genes.items():
                    data.loc[k1, k2] = len(set(v1) & set(v2))
            self._interaction_intersection = data
        return self._interaction_intersection

    @property
    def interaction_union(self) -> pd.DataFrame:
        if self._interaction_union is None:
            data = pd.DataFrame(
                columns=list(self.accession_genes.keys()),
                index=list(self.accession_genes.keys()),
                dtype=int,
            )
            for (k1, v1) in tqdm(self.accession_genes.items()):
                for (k2, v2) in self.accession_genes.items():
                    data.loc[k1, k2] = len(set(v1) | set(v2))
            self._interaction_union = data
        return self._interaction_union

    @property
    def rbp_interaction_count(self) -> pd.Series:
        return pd.Series({k: len(v) for (k, v) in self.accession_genes.items()})

    @property
    def nunique_gene(self) -> int:
        if self._nunique_gene is None:
            self._nunique_gene = len(set(chain(*self.accession_genes.values())))
        return self._nunique_gene


class Default(SimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        return super().execute()
