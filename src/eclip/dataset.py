from typing import List, Union

import pandas as pd

from src.util.bed_format_strategy import FormatStrategy
from src.util.bedfile import load_report, read_annotated_bed
from src.util.get_bed_path import get_formatted_file_path, get_file_path
from src.eclip.uniprot.keyword import Keyword


class Dataset:
    def __init__(self, dataset: Union[str, pd.Series], IDR_replicate: str = "1,2"):
        if isinstance(dataset, str):
            self._dataset: pd.Series = (
                load_report()
                .query(
                    "Dataset == @dataset & `Biological replicates` == @IDR_replicate"
                )
                .iloc[0]
            )
        else:
            self._dataset = dataset

        self._keywords: List[str]
        self._genes: List[str]

    @property
    def protein(self):
        return self._dataset["Target label"]

    @property
    def biosample(self):
        return self._dataset["Biosample name"]

    @property
    def dataset(self):
        return self._dataset["Dataset"]

    @property
    def keywords(self):
        if hasattr(self, "_keywords"):
            return self._keywords

        return Keyword()(self.protein)

    @property
    def genes(self):
        if hasattr(self, "_genes"):
            return self._genes
        df = read_annotated_bed(
            get_formatted_file_path(self._dataset, FormatStrategy.MAX)
        )
        self._genes = list(set(df["gene_name"]))
        return self._genes

    @property
    def peaks_len(self):
        peakfile = get_file_path(self._dataset)
        with open(peakfile, "r") as f:
            return len(f.readlines())

    def __repr__(self) -> str:
        return f"{self.dataset} ({self.protein}, {self.biosample})"
