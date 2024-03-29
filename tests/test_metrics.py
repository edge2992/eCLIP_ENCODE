from src.util.metrics.metrics import Metrics


def test_metrics_tape():
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
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_strategy import TAPE, Jaccard, Dice
    import pandas as pd

    TEST_N = 10

    report = load_replicateIDR_report().head(TEST_N)
    metrics = Metrics(report)
    data_list = []
    for strategy in [TAPE(symmetric_method=None), Jaccard(), Dice()]:
        data_list.append(metrics(strategy, add_description=False))
    data = pd.concat([metrics.description(), pd.DataFrame(data_list).T], axis=1)

    multiple_metrics = metrics(
        [TAPE(symmetric_method=None), Jaccard(), Dice()], add_description=True
    )

    assert all(data == multiple_metrics)

    expected_tape = metrics(TAPE(), add_description=True)
    expected_dice = metrics(Dice(), add_description=True)
    for index, row in data.iterrows():
        assert row["TAPE Cosine"] == expected_tape.loc[index, "TAPE Cosine"]  # type: ignore
        assert row["Gene Dice"] == expected_dice.loc[index, "Gene Dice"]  # type: ignore


def test_metrics_peak(sample_report):
    from src.util.similarity_strategy import PeakStrategy

    metrics = Metrics(sample_report)
    data_peak = metrics(PeakStrategy(), add_description=True)
    print(data_peak)


def test_metrics_peak_jaccard(K562_report):
    from src.util.similarity_strategy import PeakStrategy, Jaccard

    metrics = Metrics(K562_report)
    data = metrics([PeakStrategy(metrics="jaccard"), Jaccard()], add_description=True)
    expected_shape = K562_report.shape[0] * (K562_report.shape[0] - 1) / 2
    assert data.shape[0] == expected_shape
    assert any(data["Gene Jaccard"].isna()) == False
    assert any(data["Peak Jaccard"].isna()) == False


def test_qcut_stringdb(sample_report):
    from src.util.similarity_strategy import DirectStringScore

    metrics = Metrics(sample_report)
    data = metrics(DirectStringScore(metrics="score", qcut=True), add_description=False)
    assert data.shape[0] == sample_report.shape[0] * (sample_report.shape[0] - 1) / 2
    v_count = data.value_counts()
    assert v_count["l"] == 38
    assert v_count["m"] == 4
    assert v_count["h"] == 3
