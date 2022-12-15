# ある条件で実験ペアを抽出した時にキーワードがエンリッチメントされる度合いをフィッシャーの正確検定で測る
import os
from typing import Dict

import pandas as pd
import scipy.stats as st
from dotenv import load_dotenv

from src.eclip import Compare, Dataset

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
LABEL_KEYWORD_REPORT_PATH = os.path.join(PROJECT_PATH, "data/label_keyword_report.csv")
CACHE_DIR = os.path.join(PROJECT_PATH, ".cache")

if not os.path.exists(CACHE_DIR):
    os.mkdir(CACHE_DIR)


class Condition:
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


class KeywordConfidence:
    def __init__(
        self,
        report: pd.DataFrame,
        label_keyword_report_path: str = LABEL_KEYWORD_REPORT_PATH,
    ):
        self.report = report
        self.label_keyword_report_path = label_keyword_report_path
        self.__keyword_report: pd.DataFrame
        self.__keyword_dict: Dict

    def keyword_confidence(
        self, keyword: str, condition: Condition, alternatives: str = "greater"
    ):
        """conditionに合う条件で実験ペアを抽出した時にkeywordがエンリッチメントされている度合いをFisherの正確検定で検定する"""
        crosstab = self.cross_tab(keyword, condition)
        return st.fisher_exact(crosstab, alternative=alternatives)

    def cross_tab(self, keyword: str, condition: Condition):
        targets = self.report[condition(self.report)]
        assert isinstance(targets, pd.DataFrame)
        target_labels = self.labels(targets)
        data = pd.DataFrame(
            {
                str(condition): self.keyword_report["label"].isin(target_labels),  # type: ignore
                keyword: self.keyword_report["keyword"] == keyword,
            }
        )
        return (
            pd.crosstab(data[str(condition)], data[keyword])
            .sort_index(axis=0, ascending=False)
            .sort_index(axis=1, ascending=False)
        )

    @property
    def keywords_dict(self) -> Dict:
        # keyword -> List[label of pair]のDictを返す
        if not hasattr(self, "__keyword_dict"):
            self.__keyword_dict = dict(
                self.keyword_report.groupby("keyword")["label"].apply(list)
            )
        return self.__keyword_dict

    @property
    def keyword_report(self) -> pd.DataFrame:
        if not hasattr(self, "__keyword_report"):
            self.__keyword_report = self.__load_label_keyword_report()
        return self.__keyword_report

    @property
    def keywords(self):
        return self.keyword_report["keyword"].unique().tolist()

    def __load_label_keyword_report(self):
        if not os.path.exists(self.label_keyword_report_path):
            self.__make_label_keyword_report().to_csv(
                self.label_keyword_report_path, index=False
            )
        cache_path = os.path.join(CACHE_DIR, f"label_keyword_{self.__hash()}")
        if not os.path.exists(cache_path):
            df = pd.read_csv(
                self.label_keyword_report_path,
            )
            target_labels = self.labels(self.report)  # noqa: F841
            data = df.query("label in @target_labels")
            data.to_csv(cache_path, index=False)
        else:
            print(f"load cache {cache_path}")

        return pd.read_csv(cache_path)

    def __make_label_keyword_report(self) -> pd.DataFrame:
        from src.util.bedfile import load_replicateIDR_report
        from src.util.metrics import Metrics

        def keyword_intersection(row: pd.Series) -> Dict:
            compare = Compare([Dataset(row["Dataset_1"]), Dataset(row["Dataset_2"])])
            return {"label": compare.label, "keyword": compare.keyword_intersection}

        print("make label_keyword_report")
        report_all = Metrics(load_replicateIDR_report()).description()
        data = report_all.apply(keyword_intersection, axis=1, result_type="expand").explode("keyword").reset_index(drop=True)  # type: ignore
        print("done")
        return data

    @classmethod
    def labels(cls, report: pd.DataFrame):
        return report.apply(
            lambda row: Compare(
                [Dataset(row["Dataset_1"]), Dataset(row["Dataset_2"])]
            ).label,
            axis=1,
        ).to_list()  # type: ignore

    def __hash(self):
        import hashlib

        query = str(self.report.iloc[:, 0].to_list())
        return hashlib.md5(query.encode()).hexdigest()


if __name__ == "__main__":
    from src.util.metrics import Metrics

    def __HepG2_1e3_report() -> pd.DataFrame:
        from src.eclip.dataset import Dataset
        from src.util.bedfile import load_replicateIDR_report

        report = load_replicateIDR_report()
        report = report[report["Biosample name"] == "HepG2"]
        report = report[
            report.apply(lambda row: len(Dataset(row).genes) >= 1e3, axis=1)
        ]
        return report.reset_index(drop=True)

    report = Metrics(__HepG2_1e3_report()).description()
    confidence = KeywordConfidence(report)
