import os
from abc import ABC, abstractmethod
from itertools import chain
from typing import Callable, Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from dotenv import load_dotenv
from tqdm import tqdm

from src.plot.util.process_report import label_protein_biosample

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]


class SimilarityStrategy(ABC):
    def __init__(self, report: Union[None, pd.DataFrame] = None):
        if report is None:
            from src.util.bedfile import load_replicateIDR_report

            self.report = load_replicateIDR_report()
        else:
            self.report = report

    @abstractmethod
    def execute(self) -> pd.DataFrame:
        raise NotImplementedError()


class SequenceSimilarityStrategy(SimilarityStrategy):

    UNIPROT_TABLE = os.path.join(PROJECT_PATH, "data/uniprot", "reviewed.tsv")

    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
    ):
        super().__init__(report)
        self.loadfile = loadfile

    def execute(self) -> pd.DataFrame:
        similarities = self._protein_similarity()
        mapping_dict = self._idmapping()

        return similarities.loc[
            list(mapping_dict.keys()), list(mapping_dict.keys())
        ].rename(index=mapping_dict, columns=mapping_dict)

    @abstractmethod
    def _protein_similarity(self) -> pd.DataFrame:
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
        label_method: Callable[[pd.DataFrame], pd.Series] = label_protein_biosample,
    ):
        super().__init__(report)
        self.loadfile = loadfile
        self.label_method = label_method

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


class MSA(SequenceSimilarityStrategy):
    def _protein_similarity(self) -> pd.DataFrame:
        similarities, labels = self._load()
        return pd.DataFrame(similarities, index=labels, columns=labels, dtype=float)

    def _load(self) -> Tuple[List[List[str]], List[str]]:
        similarities: List[List[str]] = []
        labels: List[str] = []

        if self.loadfile is None:
            self.loadfile = os.environ["CLUSTALO_DISTMAT_PATH"]

        if not os.path.exists(self.loadfile):
            self._execute_msa(self.loadfile)

        with open(self.loadfile, "r") as f:
            lines = f.readlines()
            for line in lines[1:]:
                if line.startswith("#"):
                    continue
                elif line == "\n":
                    continue
                else:
                    similarities.append(line.split()[1:])
                    labels.append(line.split()[0])
        return similarities, labels

    def _execute_msa(self, path: str):
        """Clustal Omegaを実行してファイルをpathに保存する"""
        raise NotImplementedError()


class TAPE(SequenceSimilarityStrategy):
    def _protein_similarity(self) -> pd.DataFrame:
        from scipy.spatial.distance import cosine, pdist, squareform

        from src.util.tape import load_embedding

        if self.loadfile is None:
            self.loadfile = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_tape.npz"

        if not os.path.exists(self.loadfile):
            self._execute_tape(self.loadfile)

        embedding = load_embedding(how="avg", filepath=self.loadfile)

        return pd.DataFrame(
            squareform(pdist(np.array(list(embedding.values()), dtype=float), cosine)),
            index=list(embedding.keys()),
            columns=list(embedding.keys()),
        )

    def _execute_tape(self, path: str):
        """TAPEを実行してファイルをpathに保存する"""
        raise NotImplementedError()


class Default(SimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        return super().execute()


class Lift(InteractionSimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        """タンパク質のリフト値を計算する。
        リフト値はマーケット分析で使われる指標
        eCLIPはタンパク質を中心に解析をするが、
        視点を変えてRNAごとに結合するタンパク質が出現する頻度を考えて、
        タンパク質同士の出現の相関を求める。
        出現の仕方が互いに一切関係なく、独立な場合にはリフト値は1となる。
        正の相関だと1より大きくなり負の相関だと1より小さくなる。
        """
        rbp_count = self.rbp_interaction_count.to_numpy()

        value = (
            self.interaction_intersection.to_numpy()
            / (rbp_count * rbp_count.reshape((-1, 1)))
            * self.nunique_gene
        )

        return pd.DataFrame(
            value,
            index=self.interaction_intersection.index,
            columns=self.interaction_intersection.columns,
            dtype=float,
        )


class Dice(InteractionSimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        rbp_count = self.rbp_interaction_count.to_numpy()

        value = (
            self.interaction_intersection.to_numpy()
            * 2
            / (rbp_count + rbp_count.reshape((-1, 1)))
        )

        return pd.DataFrame(
            value,
            index=self.interaction_intersection.index,
            columns=self.interaction_intersection.columns,
            dtype=float,
        )


class Jaccard(InteractionSimilarityStrategy):
    def execute(self) -> pd.DataFrame:

        value = (
            self.interaction_intersection.to_numpy() / self.interaction_union.to_numpy()
        )

        return pd.DataFrame(
            value,
            index=self.interaction_intersection.index,
            columns=self.interaction_intersection.columns,
            dtype=float,
        )


class Simpson(InteractionSimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        rbp_count = self.rbp_interaction_count.to_numpy()
        lower = np.minimum(rbp_count, rbp_count.reshape((-1, 1)))
        assert self.interaction_intersection.shape == lower.shape

        value = self.interaction_intersection.to_numpy() / lower

        return pd.DataFrame(
            value,
            index=self.interaction_intersection.index,
            columns=self.interaction_intersection.columns,
            dtype=float,
        )


class Cosine(InteractionSimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        rbp_count = self.rbp_interaction_count.to_numpy()

        value = self.interaction_intersection.to_numpy() / np.sqrt(
            rbp_count * rbp_count.reshape((-1, 1))
        )

        return pd.DataFrame(
            value,
            index=self.interaction_intersection.index,
            columns=self.interaction_intersection.columns,
            dtype=float,
        )