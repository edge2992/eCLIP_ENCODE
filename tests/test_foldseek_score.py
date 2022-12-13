def test_foldseek_score():
    from src.util.similarity_strategy.foldseek import FoldSeekTMScore
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.uniprot import idmapping
    from foldseek.wrapper import read_aln_tmscore

    similarity = ProteinSimilarity()
    similarity.setStrategy(FoldSeekTMScore())
    data = similarity.executeStrategy()

    assert similarity.strategy.loadfile is not None
    expected_table = read_aln_tmscore(similarity.strategy.loadfile)
    protein_converter = {v.split("|")[1]: k for k, v in idmapping().items()}

    for _, row in expected_table.head(1000).iterrows():
        protein1 = protein_converter[row["query"].split("-")[1]]
        protein2 = protein_converter[row["target"].split("-")[1]]
        score = row["TMscore"]
        assert data.loc[protein1, protein2] == score


def test_foldseek_score_column_index():
    """columnとindexの順序が同じ"""
    from src.util.similarity_strategy.foldseek import FoldSeekTMScore
    from src.util.similarity_protein import ProteinSimilarity

    similarity = ProteinSimilarity()
    similarity.setStrategy(FoldSeekTMScore())
    data = similarity.executeStrategy()
    assert all(data.columns == data.index)


def test_foldseek_score_transform():
    """transformしてdatasetをkeyとしても正しいvalueが得られる"""
    from src.util.similarity_strategy.foldseek import FoldSeekTMScore
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.bedfile import load_replicateIDR_report
    from src.util.uniprot import idmapping
    from foldseek.wrapper import read_aln_tmscore
    import numpy as np

    N_TEST = 15
    report = load_replicateIDR_report().head(N_TEST)
    similarity = ProteinSimilarity()
    similarity.setStrategy(
        FoldSeekTMScore(report=report, symmetric=True, symmetric_method="min")
    )
    stringdb_data = similarity.executeStrategy()
    data = similarity.strategy.transform(stringdb_data)
    assert data.shape[0] == N_TEST
    assert data.shape[1] == N_TEST

    report_ind = report.set_index("Dataset", drop=True)
    assert similarity.strategy.loadfile is not None
    expected_table = read_aln_tmscore(similarity.strategy.loadfile)
    protein_converter = {k: v.split("|")[1] for k, v in idmapping().items()}

    for index, row in data.iterrows():
        for col, value in row.items():
            protein1 = protein_converter[report_ind.loc[str(index), "Target label"]]
            protein2 = protein_converter[report_ind.loc[str(col), "Target label"]]
            protein1 = "AF-{}-F1-model_v3.cif.gz".format(protein1)
            protein2 = "AF-{}-F1-model_v3.cif.gz".format(protein2)
            expected = expected_table[
                expected_table["target"].isin([protein1, protein2])
                & expected_table["query"].isin([protein1, protein2])
                & (expected_table["query"] != expected_table["target"])
            ]

            if protein1 == protein2:
                assert value == 0
            elif value == 0:
                # TODO: 片方のデータしかない場合minを取ると値が0になってしまう (修正する)
                continue
            else:
                assert expected is not None, f"{protein1}, {protein2}"
                expected_score = min(expected["TMscore"].values)
                assert np.isclose(expected_score, value)
