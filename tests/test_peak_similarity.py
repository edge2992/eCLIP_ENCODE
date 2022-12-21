def test_peak_similarity():
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_strategy import Simpson, PeakStrategy
    from src.util.similarity_protein import Similarity, InteractionSimilarity

    N_TEST = 20
    report = load_replicateIDR_report().head(N_TEST)
    similarity = Similarity()
    similarity.setStrategy(PeakStrategy(report=report))
    peak_data = similarity.executeStrategy()
    assert peak_data.shape == (N_TEST, N_TEST)

    similarity = InteractionSimilarity()
    similarity.setStrategy(Simpson(report=report))
    simpson_data = similarity.executeStrategy()
    assert all(peak_data.columns == simpson_data.columns)
    assert all(peak_data.index == simpson_data.index)
