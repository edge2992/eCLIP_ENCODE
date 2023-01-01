def test_keyword(sample_report):
    from src.eclip.uniprot.keyword import Keyword
    from src.eclip import Dataset

    # FMR1 35
    # SAFB 14
    # IGF2BP3 13

    keyword_search = Keyword()

    assert len(keyword_search("FMR1")) == 35
    assert len(keyword_search("SAFB")) == 14
    assert len(keyword_search("IGF2BP3")) == 13

    for _, row in sample_report.iterrows():
        keys = keyword_search(row["Target label"])
        expected = Dataset(row["Dataset"]).keywords
        assert len(keys) == len(expected)


def test_protein_keyword_table():
    from src.eclip.uniprot.keyword import Keyword

    df = Keyword.protein_keyword_table()
    assert len(df[df.Protein == "FMR1"]) == 35
    assert len(df[df.Protein == "SAFB"]) == 14
    assert len(df[df.Protein == "IGF2BP3"]) == 13

    assert len(df[df.Keyword == "3D-structure"]) == 123
    assert len(df[df.Keyword == "mRNA processing"]) == 52
    assert len(df[df.Keyword == "Transcription"]) == 38
