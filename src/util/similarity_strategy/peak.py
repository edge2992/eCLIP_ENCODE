from typing import Callable, Union

import pandas as pd
import os
from dotenv import load_dotenv

from src.util.similarity_strategy.interface import SimilarityStrategy
from src.bed_jaccard import sort_all, jaccard_bed_all, convert_df_from_jaccard_output
from src.util.bedfile import load_replicateIDR_report
from src.util.get_bed_path import get_file_path

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

CACHE_DIR = os.path.join(PROJECT_PATH, ".cache")
SORTED_DIR = os.path.join(PROJECT_PATH, "data/bed/sorted")
JACCARD_DIR = os.path.join(PROJECT_PATH, "data/bed/jaccard")

if not os.path.exists(CACHE_DIR):
    os.makedirs(CACHE_DIR)


class PeakStrategy(SimilarityStrategy):
    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
        label_method: Callable[[pd.DataFrame], pd.Series] = lambda df: df["Dataset"],
        metrics: str = "jaccard",
    ):
        super().__init__(report, loadfile, label_method)
        self.metrics = metrics

    def execute(self) -> pd.DataFrame:
        data = self._load()
        mapping_dict = self._idmapping()
        return data.loc[list(mapping_dict.keys()), list(mapping_dict.keys())].rename(
            index=mapping_dict, columns=mapping_dict
        )

    def _idmapping(self):
        df = pd.DataFrame(
            {
                "From": self.report["Accession"],
                "To": self.label_method(self.report),
            }
        )
        return df.set_index("From").to_dict()["To"]

    def _load(self) -> pd.DataFrame:
        cache_path = os.path.join(CACHE_DIR, f"peak_similarity_{self.metrics}.csv")
        if not os.path.exists(cache_path):
            self._run_bedtools()
            df = convert_df_from_jaccard_output(JACCARD_DIR, self.metrics)
            print("Saving cahce...", cache_path)
            df.to_csv(cache_path)
        return pd.read_csv(cache_path, index_col=0)

    @classmethod
    def _run_bedtools(cls):

        report = load_replicateIDR_report().sort_values("Dataset")
        bedfiles = report.apply(lambda row: get_file_path(row), axis=1).tolist()
        sorted_bedfiles = [
            os.path.join(SORTED_DIR, os.path.basename(bedfile)) for bedfile in bedfiles
        ]

        if not os.path.exists(SORTED_DIR):
            os.makedirs(SORTED_DIR)

        sort_all(bedfiles, sorted_bedfiles)

        if not os.path.exists(JACCARD_DIR):
            os.makedirs(JACCARD_DIR)

        jaccard_bed_all(sorted_bedfiles, JACCARD_DIR)
