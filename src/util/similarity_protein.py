from src.util.similarity_strategy import Default, SimilarityStrategy


class Similarity:
    stragety: SimilarityStrategy = Default()

    def setStrategy(self, strategy: SimilarityStrategy):
        if strategy is not None:
            self.stragety = strategy
        else:
            self.stragety = Default()

    def executeStrategy(self):
        return self.stragety.execute()
