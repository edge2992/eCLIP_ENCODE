def test_prepare_keywords_dict(sample_metrics):
    from src.util.metrics.keyword import KeywordConfidence

    # 3D-structure 36
    # Acetylation 15
    # Activator 1
    # Alternative splicing 28
    confidence = KeywordConfidence(sample_metrics)
    keywords_dict = confidence.keywords_dict
    assert len(keywords_dict["3D-structure"]) == 36
    assert len(keywords_dict["Acetylation"]) == 15
    assert len(keywords_dict["Activator"]) == 1
    assert len(keywords_dict["Alternative splicing"]) == 28


def test_contigency_table(HepG2_1e3_metrics):
    from src.util.metrics import KeywordConfidence, ConditionGt
    import numpy as np

    confidence = KeywordConfidence(HepG2_1e3_metrics)
    crosstab = confidence.cross_tab("Acetylation", ConditionGt("Gene Jaccard", 0.3))
    assert (crosstab.values == np.array([[25, 45], [410, 601]])).all()
    #      0    1
    # 0   25   45
    # 1  410  601


def test_keyword_confidence(HepG2_1e3_metrics):
    from src.util.metrics import KeywordConfidence, ConditionGt

    confidence = KeywordConfidence(HepG2_1e3_metrics)
    oddsratio, p_value = confidence.keyword_confidence(
        "Acetylation", ConditionGt("Gene Jaccard", 0.3)
    )
    print(p_value)
    oddsratio, p_value = confidence.keyword_confidence(
        "Acetylation", ConditionGt("Gene Jaccard", 0.3), "two-sided"
    )
    print(p_value)
