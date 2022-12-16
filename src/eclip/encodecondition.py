import pandas as pd
from src.util.metrics.condition import Condition, ConditionAnd, ConditionEq
from src.eclip.dataset import Dataset
from abc import ABC


class ENCODECondition(Condition, ABC):
    """ENCODEの実験の抽出条件"""

    pass


class GeneNumberCondition(ENCODECondition, ABC):
    def __init__(self, threshold):
        self.threshold = threshold

    def gene_number(self, data: pd.DataFrame) -> pd.Series:
        return data.apply(lambda row: len(Dataset(row).genes), axis=1)


class GeneNumberLtCondition(GeneNumberCondition):
    def __call__(self, data: pd.DataFrame) -> pd.Series:
        return self.gene_number(data) < self.threshold

    def __repr__(self) -> str:
        return f"Gene Number under {self.threshold}"


class GeneNumberGteCondition(GeneNumberCondition):
    def __call__(self, data: pd.DataFrame) -> pd.Series:
        return self.gene_number(data) >= self.threshold

    def __repr__(self) -> str:
        return f"Gene Number over {self.threshold}"


ECLIP_SAMPLESETS = [
    ConditionAnd(
        [ConditionEq("Biosample name", "HepG2"), GeneNumberGteCondition(1000)]
    ),
    ConditionAnd([ConditionEq("Biosample name", "HepG2"), GeneNumberLtCondition(1000)]),
    ConditionAnd([ConditionEq("Biosample name", "K562"), GeneNumberGteCondition(1000)]),
    ConditionAnd([ConditionEq("Biosample name", "K562"), GeneNumberLtCondition(1000)]),
]
