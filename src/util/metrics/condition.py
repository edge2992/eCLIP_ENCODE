from abc import ABC
import pandas as pd
from typing import List


class Condition(ABC):
    def __init__(self, hue, threshold):
        self.hue = hue
        self.threshold = threshold

    def __call__(self, row):
        raise NotImplementedError()


class ConditionGt(Condition):
    def __init__(self, hue: str, threshold: float):
        super().__init__(hue, threshold)

    def __call__(self, data: pd.DataFrame) -> pd.Series:
        return data[self.hue] > self.threshold

    def __repr__(self):
        return f"{self.hue} > {self.threshold}"


class ConditionLt(Condition):
    def __init__(self, hue: str, threshold: float):
        super().__init__(hue, threshold)

    def __call__(self, data: pd.DataFrame) -> pd.Series:
        return data[self.hue] < self.threshold

    def __repr__(self):
        return f"{self.hue} < {self.threshold}"


class ConditionAnd(Condition):
    def __init__(self, conditions: List[Condition]):
        self.conditions = conditions

    def __call__(self, data: pd.DataFrame) -> pd.Series:
        return pd.concat(
            [condition(data) for condition in self.conditions], axis=1
        ).all(axis=1)

    def __repr__(self):
        return f"{' and '.join([str(condition) for condition in self.conditions])}"
