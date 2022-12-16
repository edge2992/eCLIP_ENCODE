import os

import numpy as np
import pandas as pd
from dotenv import load_dotenv

from src.util.similarity_strategy import InteractionSimilarityStrategy

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]


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

    def __repr__(self) -> str:
        return "Lift"


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

    def __repr__(self) -> str:
        return "Dice"


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

    def __repr__(self) -> str:
        return "Jaccard"


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

    def __repr__(self) -> str:
        return "Simpson"


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

    def __repr__(self) -> str:
        return "Cosine"
