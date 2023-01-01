from typing import Callable, Union
import pandas as pd
from foldseek.wrapper import read_aln_tmscore
from src.util.similarity_strategy.interface import ProteinSimilarityStrategy
import os
from dotenv import load_dotenv

from src.util.uniprot import idmapping

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]

TARGET_DB_PATH = os.path.join(PROJECT_PATH, "data/afdb/rdb")
QUERY_DB_PATH = os.path.join(PROJECT_PATH, "data/afdb/rdb")


def idmapping_mmcif():
    uniprot_dict = idmapping()
    result = {}
    for key, value in uniprot_dict.items():
        result[key] = "AF-{}-F1-model_v3.cif.gz".format(value.split("|")[1])
    return result


class FoldSeekTMScore(ProteinSimilarityStrategy):
    """FoldSeekTMScore
    値が存在しない場合はTMscoreを0にする
    """

    def __init__(
        self,
        report: Union[None, pd.DataFrame] = None,
        loadfile: Union[None, str] = None,
        symmetric: bool = False,
        symmetric_method: str = "avg",
        label_method: Callable[[pd.DataFrame], pd.Series] = lambda df: df["Dataset"],
        e_value: float = 1e9,
    ):
        super().__init__(report, loadfile, symmetric, symmetric_method, label_method)
        self.e_value = e_value
        if self.loadfile is None:  # type: ignore
            self.loadfile = os.path.join(
                PROJECT_PATH, "data/afdb/aln_tmscore_{}.tsv".format(self.e_value)
            )

    def _load(self) -> pd.DataFrame:
        if not os.path.exists(self.loadfile):
            # 値を計算する
            from foldseek.wrapper import FoldSeek

            fs = FoldSeek(TARGET_DB_PATH, QUERY_DB_PATH, self.loadfile)
            fs.exec(self.e_value, create_db=True)

        return read_aln_tmscore(self.loadfile)

    def _protein_similarity(self) -> pd.DataFrame:
        data = self._load()
        matrix = data.pivot(index="query", columns="target", values="TMscore")
        assert matrix.shape[0] == matrix.shape[1]
        for i in range(matrix.shape[0]):
            matrix.iloc[i, i] = 1
        matrix.fillna(0, inplace=True)
        return matrix

    def _idmapping(self):
        """foldseek column name -> eCLIP protein name"""
        return {v: k for k, v in idmapping_mmcif().items()}

    def __repr__(self) -> str:
        return f"foldseek_tmscore_{self.symmetric_method}"

    @property
    def lower_better(self) -> bool:
        return False
