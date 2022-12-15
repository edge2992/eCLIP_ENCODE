def test_metrics_tape():
    from src.util.metrics.metrics import Metrics
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_strategy import TAPE

    TEST_N = 10

    report = load_replicateIDR_report().head(TEST_N)
    metrics = Metrics(report)
    data_array = metrics(TAPE(), add_description=False)
    assert data_array.shape[0] == TEST_N * (TEST_N - 1) / 2

    data_df = metrics(TAPE(), add_description=True)
    assert data_df.shape[0] == TEST_N * (TEST_N - 1) / 2


def test_metrics_multiple():
    from src.util.metrics.metrics import Metrics
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_strategy import TAPE, Jaccard, Dice
    import pandas as pd

    TEST_N = 10

    report = load_replicateIDR_report().head(TEST_N)
    metrics = Metrics(report)
    data_list = []
    for strategy in [TAPE(), Jaccard(), Dice()]:
        data_list.append(metrics(strategy, add_description=False))
    data = pd.concat([metrics.description(), pd.DataFrame(data_list).T], axis=1)

    multiple_metrics = metrics([TAPE(), Jaccard(), Dice()], add_description=True)

    assert all(data == multiple_metrics)

    expected_tape = metrics(TAPE(), add_description=True)
    expected_dice = metrics(Dice(), add_description=True)
    for index, row in data.iterrows():
        assert row["TAPE"] == expected_tape.loc[index, "TAPE"]  # type: ignore
        assert row["Dice"] == expected_dice.loc[index, "Dice"]  # type: ignore


def test_metrics():
    from src.util.metrics.metrics import Metrics
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_strategy import TAPE, Jaccard, Simpson
    from src.plot.interaction_metrics.representative import (
        metrics,
        similarity_strategy_dict,
    )

    TEST_N = 10

    report = load_replicateIDR_report().head(TEST_N)
    sample = metrics(report, *similarity_strategy_dict())
    data = (
        Metrics(report)([TAPE(), Jaccard(), Simpson()], add_description=True)
        .sort_values("TAPE")
        .reset_index(drop=True)
    )
    for index, row in data.iterrows():
        assert row["TAPE"] == sample.loc[index, "TAPE"]
        assert row["Simpson"] == sample.loc[index, "simpson"]