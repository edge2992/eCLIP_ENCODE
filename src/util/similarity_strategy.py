import os
from abc import ABC, abstractmethod
from typing import Dict, List, Tuple, Union

import numpy as np
import pandas as pd
from dotenv import load_dotenv


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
    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
    ):
        super().__init__(report)
        self.loadfile = loadfile

    def _load_count_data(self):
        # 相互作用数と遺伝子数をそれぞれ読み込む
        # TODO: hash化してキャッシュしたい
        # DataFrameがhash化出来ない
        # そこまで重い処理でもないので、一旦保留
        # @lru_cache(maxsize=5)
        def cached_load(
            report: pd.DataFrame,
        ) -> Tuple[np.ndarray, np.ndarray, List[str], int]:
            from src.plot.util.process_intersect_gene import count_interection
            from src.plot.util.process_report import (
                count_gene,
                gene_ids_eCLIP,
                label_protein_biosample,
            )

            label_method = label_protein_biosample

            # todo typo
            interaction_intersect: pd.DataFrame = count_interection(
                report, label_method
            )
            rbp_columns: List[str] = list(interaction_intersect.columns)
            rbp_interaction: pd.Series = count_gene(report, label_method)[rbp_columns]
            nunique_gene = len(gene_ids_eCLIP(report))
            assert len(rbp_interaction) == len(rbp_columns)
            assert (rbp_interaction.index == rbp_columns).all()
            return (
                interaction_intersect.to_numpy(),
                rbp_interaction.to_numpy(),
                rbp_columns,
                nunique_gene,
            )

        return cached_load(self.report)


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
        (
            interaction_intersect,
            rbp_interaction,
            rbp_columns,
            nunique_gene,
        ) = self._load_count_data()

        value = (
            interaction_intersect
            / (rbp_interaction * rbp_interaction.reshape((-1, 1)))
            * nunique_gene
        )

        return pd.DataFrame(value, index=rbp_columns, columns=rbp_columns, dtype=float)


class Dice(SimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        return super().execute()


class Jaccard(SimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        return super().execute()


class Simpson(SimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        return super().execute()


class Cosine(SimilarityStrategy):
    def execute(self) -> pd.DataFrame:
        return super().execute()
