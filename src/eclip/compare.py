from typing import List
from src.eclip.dataset import Dataset


class Compare:
    def __init__(self, datasets: List[Dataset]):
        assert len(datasets) == 2
        self._datasets = datasets

    @property
    def keyword_intersection(self):
        return list(set(self._datasets[0].keywords) & set(self._datasets[1].keywords))

    @property
    def gene_intersection(self):
        return list(set(self._datasets[0].genes) & set(self._datasets[1].genes))

    def __repr__(self) -> str:
        return f"{self._datasets[0]} & {self._datasets[1]}"
