def test_convert_dict_keyword():
    import pandas as pd
    from src.plot.interaction_metrics.representative import (
        target_report,
        convert_to_dict_exp_pair_by_keyword,
    )
    from src.util.metrics import Metrics
    from src.util.similarity_strategy import (
        TAPE,
        KeywordCosine,
        BlastP,
        Simpson,
        Lift,
        Cosine,
    )

    report = target_report(2500, "HepG2")
    data: pd.DataFrame = Metrics(report)(
        [
            TAPE(),
            KeywordCosine(),
            BlastP(symmetric_method="avg"),
            Simpson(),
            Lift(),
            Cosine(),
        ]
    )  # type: ignore
    keyword_dict = convert_to_dict_exp_pair_by_keyword(data)
    assert len(keyword_dict["RNA-binding"]) == 45
    assert len(keyword_dict["3D-structure"]) == 78
    assert len(keyword_dict["WD repeat"]) == 1
