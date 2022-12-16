import pandas as pd
from src.util.metrics.condition import Condition
from src.util.bedfile import load_replicateIDR_report


class SampleSetECLIP:
    def __init__(self, condition: Condition):
        self.condition = condition
        self.__reprot = self.__load(self.condition)

    def __load(self, condition: Condition) -> pd.DataFrame:
        report = load_replicateIDR_report()
        result = report[condition(report)]
        assert isinstance(result, pd.DataFrame)
        return result

    @property
    def report(self):
        return self.__reprot

    def __repr__(self) -> str:
        return f"{self.condition}".replace("=", "").replace(" ", "_")


if __name__ == "__main__":
    from src.util.metrics.condition import ConditionAnd, ConditionEq
    from src.eclip.encodecondition import GeneNumberGteCondition
    from typing import List

    conditions: List[Condition] = [
        ConditionEq("Biosample name", "HepG2"),
        GeneNumberGteCondition(1000),
    ]  # type: ignore

    sampleset = SampleSetECLIP(ConditionAnd(conditions))
    print(sampleset)
    print(sampleset.report.head())
    print(sampleset.report.shape)
