from abc import ABC
from src.util.similarity_strategy import (
    Default,
    InteractionSimilarityStrategy,
    ProteinSimilarityStrategy,
    SimilarityStrategy,
)


class Similarity(ABC):
    strategy: SimilarityStrategy = Default()

    def setStrategy(self, strategy: SimilarityStrategy):
        if strategy is not None:
            self.strategy = strategy
        else:
            self.stragety = Default()

    def executeStrategy(self):
        return self.strategy.execute()


class InteractionSimilarity(Similarity):
    strategy: InteractionSimilarityStrategy


class ProteinSimilarity(Similarity):
    strategy: ProteinSimilarityStrategy

    def executeStrategy(self, transform: bool = False):
        tmp = super().executeStrategy()
        if transform:
            return self.strategy.transform(tmp)
        return tmp
