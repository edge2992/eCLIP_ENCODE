from abc import ABC
from src.util.similarity_strategy import (
    Default,
    InteractionSimilarityStrategy,
    ProteinSimilarityStrategy,
    SimilarityStrategy,
    PeakStrategy,
)
from typing import Union
import numpy as np
import pandas as pd


class Similarity(ABC):
    strategy: SimilarityStrategy = Default()

    def setStrategy(
        self, strategy: SimilarityStrategy, report: Union[None, pd.DataFrame] = None
    ):
        if strategy is not None:
            self.strategy = strategy
        else:
            self.strategy = Default()

        if report is not None:
            self.strategy.set_report(report)

    def executeStrategy(self):
        return self.strategy.execute()

    @classmethod
    def flatten_tri(cls, X: pd.DataFrame, include_diagonal: bool = True):
        """
        三角行列 + 対角行列を一次元配列に変換する
        """
        assert X.shape[0] == X.shape[1]
        return X.to_numpy()[np.triu_indices(X.shape[0], k=1 - include_diagonal)]


class InteractionSimilarity(Similarity):
    strategy: InteractionSimilarityStrategy


class PeakSimilarity(Similarity):
    strategy: PeakStrategy


class ProteinSimilarity(Similarity):
    strategy: ProteinSimilarityStrategy

    def executeStrategy(self, transform: bool = False):
        tmp = super().executeStrategy()
        if transform:
            return self.strategy.transform(tmp)
        return tmp
