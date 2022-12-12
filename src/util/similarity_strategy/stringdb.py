from typing import Union
import pandas as pd

from src.util.similarity_strategy import ProteinSimilarityStrategy
from src.util.download.stringdb import (
    ENCODEprotein2preferredName,
    download_interaction_partners,
)


class DirectStringScore(ProteinSimilarityStrategy):
    """StringDB direct score
    同じタンパク質の場合スコアを1とする。
    コアがない場合は0とする
    """

    def __init__(
        self, report: Union[None, pd.DataFrame] = None, metrics: str = "score"
    ):
        super().__init__(report=report)
        assert metrics in [
            "score",
            "nscore",
            "fscore",
            "pscore",
            "ascore",
            "escore",
            "dscore",
            "tscore",
        ], "invalid metrics {metrics}"
        self.metrics = metrics

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
        matrix.fillna(0, inplace=True)
        return matrix

    def _idmapping(self):
        return {v: k for k, v in ENCODEprotein2preferredName().items()}