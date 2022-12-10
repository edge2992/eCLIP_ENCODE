import pandas as pd

from src.util.similarity_strategy import ProteinSimilarityStrategy
from src.util.download.stringdb import (
    ENCODEprotein2preferredName,
    download_interaction_partners,
)


class DirectStringScore(ProteinSimilarityStrategy):
    """StringDB direct score"""

    def _load(self, required_score: int = 0) -> pd.DataFrame:

        return download_interaction_partners(required_score=required_score)

    def _protein_similarity(self) -> pd.DataFrame:
        data = self._load()
        target_preferredNames = data["preferredName_A"].unique()
        data = data[data["preferredName_B"].isin(target_preferredNames)]
        return data.pivot(
            index="preferredName_A", columns="preferredName_B", values="score"
        )

    def _idmapping(self):
        return {v: k for k, v in ENCODEprotein2preferredName().items()}
