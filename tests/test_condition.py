import pandas as pd
import numpy as np
from src.util.metrics import Metrics
from src.util.similarity_strategy import DirectStringScore


def test_condition_neq(HepG2_1e3_report):
    from src.util.metrics.condition import ConditionNeq

    data = Metrics(HepG2_1e3_report)(
        DirectStringScore(symmetric_method=None), add_description=True
    )
    assert isinstance(data, pd.DataFrame)
    print(data.head())
    non_zero_data = data[ConditionNeq("STRING Score", 0)(data)]
    assert (data["STRING Score"] == 0).any()
    assert (non_zero_data["STRING Score"] != 0).all()


def test_condition_ltquantile(HepG2_1e3_nonzero_stringdb_metrics):
    from src.util.metrics.condition import ConditionLtQuantile

    sample = HepG2_1e3_nonzero_stringdb_metrics[
        ConditionLtQuantile("STRING Score", 0.25)(HepG2_1e3_nonzero_stringdb_metrics)
    ]

    assert np.isclose(
        sample.shape[0] / HepG2_1e3_nonzero_stringdb_metrics.shape[0], 0.25, rtol=0.1
    )


def test_condition_gtquantile(HepG2_1e3_nonzero_stringdb_metrics):
    from src.util.metrics.condition import ConditionGtQuantile

    sample = HepG2_1e3_nonzero_stringdb_metrics[
        ConditionGtQuantile("STRING Score", 0.75)(HepG2_1e3_nonzero_stringdb_metrics)
    ]

    assert np.isclose(
        sample.shape[0] / HepG2_1e3_nonzero_stringdb_metrics.shape[0], 0.25, rtol=0.1
    )
