def test_stringdb_direct_similarity():
    from src.util.similarity_strategy.stringdb import DirectStringScore
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.download.stringdb import (
        download_interaction_partners,
        ENCODEprotein2preferredName,
    )

    similarity = ProteinSimilarity()
    similarity.setStrategy(DirectStringScore())
    data = similarity.executeStrategy()
    print(data.head())

    original_table = download_interaction_partners()
    converter = ENCODEprotein2preferredName()
    rev_converter = {v: k for k, v in converter.items()}

    expected_table = original_table[
        original_table["preferredName_B"].isin(converter.values())
    ]

    for _, row in expected_table.head(1000).iterrows():
        protein1 = rev_converter[row["preferredName_A"]]
        protein2 = rev_converter[row["preferredName_B"]]
        assert data.loc[protein1, protein2] == row["score"]
