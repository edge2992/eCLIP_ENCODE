from src.util.metrics.protein_metrics import ProteinMetrics


def test_proteins():
    metrics = ProteinMetrics()
    proteins = metrics.proteins()
    assert len(proteins) == 159


def test_description():
    metrics = ProteinMetrics()
    df = metrics.description()
    assert df.shape == (159 * 158 / 2, 2)
    assert df["Protein_1"].nunique() == 158
    assert df["Protein_2"].nunique() == 158


def test_similarity():
    from src.util.similarity_strategy import TAPE
    from src.util.similarity_protein import ProteinSimilarity

    metrics = ProteinMetrics()
    strategy = TAPE()
    data = metrics(strategy, add_description=True)
    print(data.head())

    handler = ProteinSimilarity()
    handler.setStrategy(strategy)
    expected = handler.executeStrategy(transform=False)

    for _, row in data.head(100).iterrows():
        protein_a = row["Protein_1"]
        protein_b = row["Protein_2"]
        assert row[str(strategy)] == expected[protein_a][protein_b]
