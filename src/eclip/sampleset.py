import pandas as pd
from typing import List
from src.util.metrics.condition import Condition
from src.util.bedfile import load_replicateIDR_report


class SampleSetECLIP:
    def __init__(self, condition: Condition, drop_IDR_duplicates: bool = True):
        self.condition = condition
        self.__report = self.__load(self.condition, drop_IDR_duplicates)

    def __load(self, condition: Condition, drop_IDR_duplicates: bool) -> pd.DataFrame:
        report = load_replicateIDR_report(drop_duplicates=drop_IDR_duplicates)
        result = report[condition(report)]  # type: ignore
        assert isinstance(result, pd.DataFrame)
        return result

    @property
    def report(self):
        return self.__report

    @property
    def proteins(self) -> List:
        return list(self.__report["Target label"].unique())

    def __repr__(self) -> str:
        return f"{self.condition}".replace("=", "").replace(" ", "_")


if __name__ == "__main__":
    from src.util.metrics.condition import ConditionAnd, ConditionEq
    from src.eclip.encodecondition import GeneNumberGteCondition

    conditions: List[Condition] = [
        ConditionEq("Biosample name", "HepG2"),
        GeneNumberGteCondition(1000),
    ]  # type: ignore

    sampleset = SampleSetECLIP(ConditionAnd(conditions))
    print(sampleset)
    print(sampleset.report.head())
    print(sampleset.report.shape)
