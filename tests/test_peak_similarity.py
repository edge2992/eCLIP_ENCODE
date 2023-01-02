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


def test_peak_shape(HepG2_report):
    from src.util.similarity_protein import PeakSimilarity
    from src.util.similarity_strategy import PeakStrategy

    handler = PeakSimilarity()
    handler.setStrategy(PeakStrategy(metrics="jaccard"), HepG2_report)
    data = handler.executeStrategy()
    print(HepG2_report.shape)
    print(data.shape)
    assert data.shape[0] == HepG2_report.shape[0]
    assert data.shape[0] == data.shape[1]
    print(data.head())
    print(handler.flatten_tri(data, include_diagonal=False).shape)
