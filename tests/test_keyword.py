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
