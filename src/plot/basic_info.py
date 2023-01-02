# Dataset Metrics x Interaction Metrics (Jaccard, Simpson) について以下の三つのグラフをそれぞれ作成する
# 1. 上位 0.25 Quantileに対するクロス集計図のヒートマップ
# 2. hueを上位0.25 Dataset Metricsとした時のInteraction Metricsのバイオリンプロット
# 3. 0.25, 0.5, 0.75 Quantileに対するFisherの正確検定の結果をヒートマップにまとめる


# 入力条件
# HepG2 Over 1000 interaction
# HepG2 Under 1000 interaction
# K562 Over 1000 interaction
# K562 Under 1000 interaction
# %%

import os
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv
from scipy.stats import fisher_exact
from tqdm import tqdm

from src.eclip.encodecondition import ECLIP_SAMPLESETS
from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import (
    Condition,
    ConditionAnd,
    ConditionGtQuantile,
    ConditionLtQuantile,
)
from src.util.similarity_strategy import (
    TAPE,
    BlastP,
    DirectStringScore,
    FoldSeekTMScore,
    Jaccard,
    KeywordCosine,
)
from src.util.similarity_strategy.interface import SimilarityStrategy

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(PROJECT_PATH, "src/plot/img", "basic_info")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# %%


class BasicInfo:
    THRESHOLDS = np.arange(0.2, 1, 0.2, dtype=float)
    ALTERNATIVE = ["two-sided", "less", "greater"]

    def __init__(
        self,
        eclip_sampleset: ConditionAnd = ECLIP_SAMPLESETS[0],
        similarity_strategies: List[SimilarityStrategy] = [
            DirectStringScore(),
            Jaccard(),
        ],
        condition_interaction: Condition = ConditionGtQuantile("Jaccard", 0.75),
        condition_metrics: Condition = ConditionGtQuantile("stringdb_score", 0.75),
    ):
        self.sampleset = SampleSetECLIP(eclip_sampleset)
        data = Metrics(self.sampleset.report)(similarity_strategies)
        assert isinstance(data, pd.DataFrame)
        self.data: pd.DataFrame = data
        self.condition_interaction = condition_interaction
        self.condition_metrics = condition_metrics

    def plot_basic(self, thresholds=np.arange(0.2, 1, 0.2, dtype=float)):
        fig, axes = plt.subplots(len(thresholds), 2, figsize=(10, 5 * len(thresholds)))
        for i, the in enumerate(thresholds):
            self.condition_metrics.set_threshold(the)
            ctab = self.crosstab(
                self.data, self.condition_metrics, self.condition_interaction
            )
            self.plot_crosstab(ctab, f"{self.condition_metrics}", ax=axes[i, 0])  # type: ignore
            self.plot_violin(self.data, self.condition_metrics, ax=axes[i, 1])  # type: ignore
        return fig

    def plot_fisher_exact(self, thresholds=np.arange(0.2, 1, 0.2, dtype=float)):
        fisher_exact_result: Dict[str, Dict] = {}
        for i, the in enumerate(thresholds):
            self.condition_metrics.set_threshold(the)
            ctab = self.crosstab(
                self.data, self.condition_metrics, self.condition_interaction
            )
            fisher_exact_result[str(self.condition_metrics)] = {}
            for alt in self.ALTERNATIVE:
                fisher_exact_result[str(self.condition_metrics)][alt] = fisher_exact(
                    ctab, alternative=alt
                )[1]

        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        sns.heatmap(pd.DataFrame(fisher_exact_result).T, annot=True, fmt=".2e", ax=ax)
        fig.suptitle("Fisher Exact Test\n" + str(self.sampleset))
        fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.95])  # type: ignore
        return fig

    @classmethod
    def qcut_condition(cls, data: pd.DataFrame, condition: Condition):
        return condition(data).replace({True: "Strong", False: "Weak"})

    @classmethod
    def crosstab(
        cls, data: pd.DataFrame, row_condition: Condition, col_condition: Condition
    ):
        return pd.crosstab(
            BasicInfo.qcut_condition(data, row_condition),
            BasicInfo.qcut_condition(data, col_condition),
        )

    @classmethod
    def plot_crosstab(cls, crosstab: pd.DataFrame, title: str, ax=None):
        return sns.heatmap(crosstab, annot=True, fmt="d", ax=ax).set_title(title)

    @classmethod
    def plot_violin(cls, data: pd.DataFrame, condition: Condition, ax=None):
        return sns.violinplot(
            data, x=BasicInfo.qcut_condition(data, condition), y="Jaccard", ax=ax
        ).set_title(f"Strong (#{condition(data).sum()}, {condition})")


# %%
# metrics, condition(GtQuantile or LtQuantile)
def input_basicinfo(strategy: SimilarityStrategy):
    condition_metrics = (
        ConditionLtQuantile(str(strategy), 0.0)
        if strategy.lower_better
        else ConditionGtQuantile(str(strategy), 0.0)
    )
    return [Jaccard(), strategy], condition_metrics


STRATEGIES = [
    TAPE(symmetric_method="avg"),
    BlastP(symmetric_method="avg"),
    KeywordCosine(),
    DirectStringScore(metrics="score"),
    DirectStringScore(metrics="ascore"),
    DirectStringScore(metrics="escore"),
    DirectStringScore(metrics="tscore"),
    FoldSeekTMScore(symmetric_method="avg"),
]


for eclip_sampleset in ECLIP_SAMPLESETS:
    for strategy in tqdm(STRATEGIES):
        strategy_list, condition_metrics = input_basicinfo(strategy)
        basicinfo = BasicInfo(
            eclip_sampleset, strategy_list, condition_metrics=condition_metrics  # type: ignore
        )
        savedir = os.path.join(SAVEDIR, str(basicinfo.sampleset))
        if not os.path.exists(savedir):
            os.makedirs(savedir)
        fig = basicinfo.plot_basic()
        fig.savefig(os.path.join(savedir, f"basicinfo_{strategy}.png"))

        fig = basicinfo.plot_fisher_exact()
        fig.savefig(os.path.join(savedir, f"fisher_exact_{strategy}.png"))
        plt.close()

# %%
