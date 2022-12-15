import pytest
import pandas as pd


@pytest.fixture
def sample_metrics() -> pd.DataFrame:
    from src.util.metrics import Metrics
    from src.util.similarity_strategy import TAPE, Jaccard
    from src.util.bedfile import load_replicateIDR_report

    TEST_N = 10
    report = load_replicateIDR_report().head(TEST_N)
    data = Metrics(report)([TAPE(), Jaccard()], add_description=True)
    assert isinstance(data, pd.DataFrame)
    return data


@pytest.fixture
def sample_report() -> pd.DataFrame:
    from src.util.bedfile import load_replicateIDR_report

    TEST_N = 10
    return load_replicateIDR_report().head(TEST_N)
