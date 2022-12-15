# ある条件で実験ペアを抽出した時にキーワードがエンリッチメントされる度合いをフィッシャーの正確検定で測る
import os
from typing import Dict

import pandas as pd
from dotenv import load_dotenv

from src.eclip import Compare, Dataset

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
LABEL_KEYWORD_REPORT_PATH = os.path.join(PROJECT_PATH, "data/label_keyword_report.csv")


class Condition:
    def __init__(self, hue, threshold):
        self.hue = hue
        self.threshold = threshold

    def __call__(self, row):
        raise NotImplementedError()


class ConditionGt(Condition):
    def __init__(self, hue: str, threshold: float):
        super().__init__(hue, threshold)

    def _execute(self, data: pd.DataFrame) -> pd.Series:
        return data[self.hue] > self.threshold

    def __repr__(self):
        return f"{self.hue} > {self.threshold}"


class ConditionLt(Condition):
    def __init__(self, hue: str, threshold: float):
        super().__init__(hue, threshold)

    def _execute(self, data: pd.DataFrame) -> pd.Series:
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

    def keyword_confidence(self, keyword: str, condition: Condition):
        """conditionに合う条件で実験ペアを抽出した時にkeywordがエンリッチメントされている度合いをp値で返す"""

        # クロス集計図を準備する

    def prepare_metrics(self, metrics: str):
        pass

    def prepare_keywords_dict(self) -> Dict:
        # keyword -> List[label of pair]のDictを返す
        df = self.__load_label_keyword_report()
        return dict(df.groupby("keyword")["label"].apply(list))

    def __load_label_keyword_report(self):
        def labels():
            return self.report.apply(
                lambda row: Compare(
                    [Dataset(row["Dataset_1"]), Dataset(row["Dataset_2"])]
                ).label,
                axis=1,
            ).to_list()

        if not os.path.exists(self.label_keyword_report_path):
            self.__make_label_keyword_report().to_csv(
                self.label_keyword_report_path, index=False
            )

        df = pd.read_csv(
            self.label_keyword_report_path,
        )
        return df.query("label in @labels()")

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


if __name__ == "__main__":
    from src.util.bedfile import load_replicateIDR_report
    from src.util.metrics import Metrics

    report = Metrics(load_replicateIDR_report()).description()
    confidence = KeywordConfidence(report)
    print(confidence.prepare_keywords_dict())
