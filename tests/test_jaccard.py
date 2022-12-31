def test_jaccard_K562(K562_report):
    from src.util.similarity_protein import InteractionSimilarity
    from src.util.similarity_strategy import Jaccard

    handler = InteractionSimilarity()
    handler.setStrategy(Jaccard(), K562_report)
    data = handler.executeStrategy()
    print(data.shape)
