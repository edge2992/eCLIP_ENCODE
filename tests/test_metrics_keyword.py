def test_prepare_keywords_dict(sample_metrics):
    from src.util.metrics.keyword import KeywordConfidence

    # 3D-structure 36
    # Acetylation 15
    # Activator 1
    # Alternative splicing 28
    confidence = KeywordConfidence(sample_metrics)
    keywords_dict = confidence.prepare_keywords_dict()
    assert len(keywords_dict["3D-structure"]) == 36
    assert len(keywords_dict["Acetylation"]) == 15
    assert len(keywords_dict["Activator"]) == 1
    assert len(keywords_dict["Alternative splicing"]) == 28
