from typing import Callable, Union
import pandas as pd

from src.util.similarity_strategy.interface import ProteinSimilarityStrategy
from src.util.download.stringdb import (
    ENCODEprotein2preferredName,
    download_interaction_partners,
)

STRINGDB_SCORE_METRICS = [
    "score",
    "escore",
    "dscore",
    "nscore",
    "fscore",
    "pscore",
    "ascore",
    "tscore",
]


class DirectStringScore(ProteinSimilarityStrategy):
    """StringDB direct score
    同じタンパク質の場合スコアを1とする。
    コアがない場合は0とする
    """

    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
        symmetric_method: Union[None, str] = None,
        label_method: Callable[[pd.DataFrame], pd.Series] = lambda df: df["Dataset"],
        fillna: Union[None, float] = 0,
        metrics: str = "score",
    ):
        super().__init__(report, loadfile, symmetric_method, label_method)
        assert metrics in STRINGDB_SCORE_METRICS, "invalid metrics {metrics}"
        self.metrics = metrics
        self.fillna = fillna

    def _load(self, required_score: int = 0) -> pd.DataFrame:

        return download_interaction_partners(required_score=required_score)

    def _protein_similarity(self) -> pd.DataFrame:
        data = self._load()
        target_preferredNames = data["preferredName_A"].unique()
        data = data[data["preferredName_B"].isin(target_preferredNames)]
        matrix = data.pivot(
            index="preferredName_A", columns="preferredName_B", values=self.metrics
        )
        assert matrix.shape[0] == matrix.shape[1]
        for i in range(matrix.shape[0]):
            matrix.iloc[i, i] = 1
        if self.fillna is not None:
            matrix.fillna(self.fillna, inplace=True)
        return matrix

    def _idmapping(self):
        return {v: k for k, v in ENCODEprotein2preferredName().items()}

    def __repr__(self) -> str:
        if self.symmetric_method is None:
            return f"STRING {self.metrics.capitalize()}"
        else:
            return f"STRING {self.metrics.capitalize()} {self.symmetric_method}"

    @property
    def lower_better(self) -> bool:
        return False
