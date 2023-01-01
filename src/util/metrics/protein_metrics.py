from typing import Union, List
import pandas as pd
import numpy as np
from src.util.similarity_strategy import ProteinSimilarityStrategy
from src.util.similarity_protein import ProteinSimilarity


class ProteinMetrics:
    def __init__(self):
        pass

    def __call__(
        self,
        strategy: Union[ProteinSimilarityStrategy, List[ProteinSimilarityStrategy]],
        add_description: bool = True,
    ) -> Union[pd.Series, pd.DataFrame]:
        if isinstance(strategy, list):
            data = pd.concat([self.__similarity(s) for s in strategy], axis=1)
        else:
            data = self.__similarity(strategy)

        if add_description:
            return pd.concat([self.description(), data], axis=1)
        else:
            return data

    def __similarity(self, strategy: ProteinSimilarityStrategy):
        data = self.__protein_similarity(strategy)
        return pd.Series(data, name=str(strategy))

    def __protein_similarity(self, strategy: ProteinSimilarityStrategy):
        handler = ProteinSimilarity()
        handler.setStrategy(strategy)
        return handler.flatten_tri(
            handler.executeStrategy(transform=False), include_diagonal=False
        )

    def description(self):
        protein_genes = self.proteins()
        N = len(protein_genes)
        index_number = np.where(np.tri(N, k=-1, dtype=bool).T)
        return pd.DataFrame(
            {
                "Protein_1": [protein_genes[i] for i in index_number[0]],
                "Protein_2": [protein_genes[i] for i in index_number[1]],
            }
        )

    @classmethod
    def proteins(cls) -> List:
        from src.util.bedfile import load_replicateIDR_report

        report = load_replicateIDR_report()
        return sorted(list(report["Target label"].unique()))
