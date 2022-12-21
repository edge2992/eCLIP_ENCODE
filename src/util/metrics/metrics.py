import numpy as np
from src.util.similarity_strategy import (
    SimilarityStrategy,
    ProteinSimilarityStrategy,
    InteractionSimilarityStrategy,
    PeakStrategy,
)
from src.util.similarity_protein import (
    ProteinSimilarity,
    InteractionSimilarity,
    PeakSimilarity,
)
import pandas as pd
from typing import Union, List


REPORT_COLUMNS = ["Dataset", "Target label", "Biosample name"]


class Metrics:
    def __init__(self, report: pd.DataFrame):
        self.report = report

    def __call__(
        self,
        strategy: Union[SimilarityStrategy, List[SimilarityStrategy]],
        add_description: bool = True,
        report_columns: list = REPORT_COLUMNS,
    ) -> Union[pd.Series, pd.DataFrame]:
        if isinstance(strategy, list):
            data = pd.concat([self.__similarity(s) for s in strategy], axis=1)
        else:
            data = self.__similarity(strategy)

        if add_description:
            return pd.concat([self.description(report_columns), data], axis=1)
        else:
            return data

    def __similarity(self, strategy: SimilarityStrategy):
        if isinstance(strategy, ProteinSimilarityStrategy):
            data = self.__protein_similarity(strategy)
        elif isinstance(strategy, InteractionSimilarityStrategy):
            data = self.__interaction_similarity(strategy)
        elif isinstance(strategy, PeakStrategy):
            data = self.__peak_similarity(strategy)
        else:
            raise ValueError(
                "strategy must be ProteinSimilarityStrategy or InteractionSimilarityStrategy"
            )
        return pd.Series(data, name=str(strategy))

    def __protein_similarity(self, strategy: ProteinSimilarityStrategy):
        handler = ProteinSimilarity()
        handler.setStrategy(strategy, self.report)
        return handler.flatten_tri(
            handler.executeStrategy(transform=True), include_diagonal=False
        )

    def __peak_similarity(self, strategy: PeakStrategy):
        handler = PeakSimilarity()
        handler.setStrategy(strategy, self.report)
        return handler.flatten_tri(handler.executeStrategy(), include_diagonal=False)

    def __interaction_similarity(self, strategy: InteractionSimilarityStrategy):
        handler = InteractionSimilarity()
        handler.setStrategy(strategy, self.report)
        return handler.flatten_tri(handler.executeStrategy(), include_diagonal=False)

    def description(self, report_columns: list = REPORT_COLUMNS):
        N = self.report.shape[0]
        index_number = np.where(np.tri(N, k=-1, dtype=bool).T)
        return pd.concat(
            [
                self.__dataset_columns(index_number[0], report_columns).add_suffix(
                    "_1"
                ),
                self.__dataset_columns(index_number[1], report_columns).add_suffix(
                    "_2"
                ),
            ],
            axis=1,
        )

    def __dataset_columns(self, index_number: np.ndarray, report_columns: list):
        return (
            self.report.iloc[index_number].loc[:, report_columns].reset_index(drop=True)
        )
