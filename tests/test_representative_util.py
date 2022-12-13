def test_convert_dict_keyword():
    from src.plot.interaction_metrics.representative import (
        metrics,
        similarity_strategy_dict,
        target_report,
        convert_to_dict_exp_pair_by_keyword,
    )

    report = target_report(2500, "HepG2")
    data = metrics(report, *similarity_strategy_dict())
    keyword_dict = convert_to_dict_exp_pair_by_keyword(data)
    assert len(keyword_dict["RNA-binding"]) == 45
    assert len(keyword_dict["3D-structure"]) == 78
    assert len(keyword_dict["WD repeat"]) == 1
