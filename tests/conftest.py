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


def __HepG2_1e3_report() -> pd.DataFrame:
    from src.util.bedfile import load_replicateIDR_report
    from src.eclip.dataset import Dataset

    report = load_replicateIDR_report()
    report = report[report["Biosample name"] == "HepG2"]
    report = report[report.apply(lambda row: len(Dataset(row).genes) >= 1e3, axis=1)]
    return report.reset_index(drop=True)


@pytest.fixture
def HepG2_1e3_report() -> pd.DataFrame:
    return __HepG2_1e3_report()


@pytest.fixture
def HepG2_1e3_metrics() -> pd.DataFrame:
    from src.util.metrics import Metrics
    from src.util.similarity_strategy import Jaccard

    data = Metrics(__HepG2_1e3_report())([Jaccard()], add_description=True)
    assert isinstance(data, pd.DataFrame)
    return data


@pytest.fixture
def HepG2_1e3_nonzero_stringdb_metrics() -> pd.DataFrame:
    from src.util.metrics import Metrics
    from src.util.metrics.condition import ConditionNeq
    from src.util.similarity_strategy import DirectStringScore

    data = Metrics(__HepG2_1e3_report())(DirectStringScore(), add_description=True)
    assert isinstance(data, pd.DataFrame)
    return data[ConditionNeq("stringdb_score", 0)(data)]


@pytest.fixture
def HepG2_report() -> pd.DataFrame:
    from src.eclip.sampleset import SampleSetECLIP
    from src.util.metrics.condition import ConditionEq

    return SampleSetECLIP(ConditionEq("Biosample name", "HepG2")).report


@pytest.fixture
def K562_report() -> pd.DataFrame:
    from src.eclip.sampleset import SampleSetECLIP
    from src.util.metrics.condition import ConditionEq

    return SampleSetECLIP(ConditionEq("Biosample name", "K562")).report
