def get_original_table():
    from src.util.download.stringdb import (
        download_interaction_partners,
        ENCODEprotein2preferredName,
    )

    original_table = download_interaction_partners(required_score=0)
    converter = ENCODEprotein2preferredName()
    expected_table = original_table[
        original_table["preferredName_B"].isin(converter.values())
    ]
    return expected_table


def test_stringdb_direct_similarity():
    from src.util.similarity_strategy.stringdb import DirectStringScore
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.download.stringdb import (
        ENCODEprotein2preferredName,
    )

    similarity = ProteinSimilarity()
    similarity.setStrategy(DirectStringScore())
    data = similarity.executeStrategy()
    print(data.head())

    expected_table = get_original_table()
    rev_converter = {v: k for k, v in ENCODEprotein2preferredName().items()}

    for _, row in expected_table.head(1000).iterrows():
        protein1 = rev_converter[row["preferredName_A"]]
        protein2 = rev_converter[row["preferredName_B"]]
        assert data.loc[protein1, protein2] == row["score"]


def test_stringdb_TAPE_column_index():
    from src.util.similarity_strategy.stringdb import DirectStringScore
    from src.util.similarity_strategy import TAPE
    from src.util.similarity_protein import ProteinSimilarity

    similarity = ProteinSimilarity()
    similarity.setStrategy(DirectStringScore())
    stringdb_data = similarity.executeStrategy()
    similarity.setStrategy(TAPE())
    tape_data = similarity.executeStrategy()
    assert all(stringdb_data.columns == tape_data.columns)


def test_stringdb_qcut():
    import pandas as pd
    from src.util.similarity_strategy.stringdb import DirectStringScore
    from src.util.similarity_protein import ProteinSimilarity

    handler = ProteinSimilarity()
    handler.setStrategy(DirectStringScore(metrics="score", qcut=True))
    data = pd.Series(
        handler.flatten_tri(handler.executeStrategy(), include_diagonal=False)
    )
    value_count = data.value_counts()
    assert value_count["low-confidence"] == 10367
    assert value_count["midium-confidence"] == 1018
    assert value_count["high-confidence"] == 1176


def test_stringdb_transform():
    from src.util.similarity_strategy.stringdb import DirectStringScore
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.bedfile import load_replicateIDR_report
    from src.util.download.stringdb import (
        ENCODEprotein2preferredName,
    )
    import numpy as np

    N_TEST = 15
    report = load_replicateIDR_report().head(N_TEST)
    similarity = ProteinSimilarity()
    similarity.setStrategy(DirectStringScore(report=report))
    stringdb_data = similarity.executeStrategy()
    data = similarity.strategy.transform(stringdb_data)
    assert data.shape[0] == N_TEST
    assert data.shape[1] == N_TEST

    report_ind = report.set_index("Dataset", drop=True)
    expected_table = get_original_table()
    converter = ENCODEprotein2preferredName()

    for index, row in data.iterrows():
        for col, value in row.items():
            protein1 = converter[report_ind.loc[str(index), "Target label"]]
            protein2 = converter[report_ind.loc[str(col), "Target label"]]
            expected = expected_table[
                expected_table["preferredName_A"].isin([protein1, protein2])
                & expected_table["preferredName_B"].isin([protein1, protein2])
            ]

            if protein1 == protein2:
                assert value == 0
            elif value == 0:
                continue
            else:
                assert not expected.empty
                expected_value = expected["score"].values[0]  # type: ignore
                assert np.isclose(expected_value, value)  # type: ignore


def test_stringdb_homology():
    from src.util.download.stringdb import (
        download_homology_metrics,
        ENCODEprotein2stringdb,
    )
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.similarity_strategy import BlastP
    import pandas as pd
    import numpy as np

    data = download_homology_metrics()
    converter = {v: k for k, v in ENCODEprotein2stringdb("stringId").items()}
    data_converted = pd.DataFrame(
        {
            "protein_1": data["stringId_A"].apply(lambda x: converter[x]),
            "protein_2": data["stringId_B"].apply(lambda x: converter[x]),
            "bitscore": data["bitscore"],
        }
    )
    print(data_converted.head())
    similarity = ProteinSimilarity()
    similarity.setStrategy(BlastP(symmetric_method="max"))
    expected = similarity.executeStrategy()
    for index, row in data_converted.iterrows():
        assert np.isclose(
            expected.loc[row["protein_1"], row["protein_2"]], row["bitscore"], 2
        )
