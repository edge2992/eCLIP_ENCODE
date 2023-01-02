import numpy as np
from src.util.similarity_protein import InteractionSimilarity
from tests.util.similarity import get_gene_list


def test_gene_n_min(sample_report):
    from src.util.similarity_strategy import Gene_N_Min

    handler = InteractionSimilarity()
    handler.setStrategy(
        Gene_N_Min(report=sample_report, label_method=lambda row: row["Accession"])
    )
    data = handler.executeStrategy()
    assert data.shape[0] == data.shape[1]
    for index, row in data.head(3).iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(sample_report, str(index))
            gene2 = get_gene_list(sample_report, str(col))
            expected = min(len(gene1), len(gene2))
            assert np.isclose(value, expected)


def test_gene_n_union(sample_report):
    from src.util.similarity_strategy import Gene_N_Union

    handler = InteractionSimilarity()
    handler.setStrategy(
        Gene_N_Union(report=sample_report, label_method=lambda row: row["Accession"])
    )
    data = handler.executeStrategy()
    assert data.shape[0] == data.shape[1]
    for index, row in data.head(3).iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(sample_report, str(index))
            gene2 = get_gene_list(sample_report, str(col))
            expected = len(set(gene1) | set(gene2))
            assert np.isclose(value, expected)
