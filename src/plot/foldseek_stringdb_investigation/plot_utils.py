from abc import ABC
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats


def construct_plotting_data(data: pd.DataFrame, keyword_experiment_pair, metrics: List):
    """plot_boxplot_by_keywordで作図するデータを作成する

    Args:
        data (pd.DataFrame): from src.plot.interaction_metrics.representativeのmetricsで作成したデータ
        keyword_experiment_pair (_type_): convert_to_dict_exp_pair_by_keywordで作成したデータ
        metrics (List): 対象とするmetrics

    Returns:
        _type_: _description_
    """
    plot_data = pd.DataFrame()
    for keyword, value in keyword_experiment_pair.items():
        sample = data.iloc[value, :][metrics].copy().reset_index()
        sample["label"] = f"{keyword} (#{len(value)})"
        plot_data = pd.concat([plot_data, sample], axis=0)
    return plot_data


def plot_boxplot_by_keyword(data: pd.DataFrame, metrics: str):
    """keywordごとにmetricsの分布をboxplotで表示する
    dataはconstruct_plotting_dataで作成したものを使用する
    """
    order_xlabels = (
        data.groupby("label")
        .mean()
        .sort_values(metrics, ascending=False)
        .index.to_list()
    )

    fig, ax = plt.subplots(figsize=(50, 10))
    sns.boxplot(data, x="label", y=metrics, order=order_xlabels, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    fig.tight_layout()
    return fig


class ConditionPlot(ABC):
    def __init__(self, hue: str, threshold: float):
        self.hue = hue
        self.threshold = threshold

    def _execute(self, data: pd.DataFrame) -> pd.Series:
        # return data[self.hue] < self.threshold
        raise NotImplementedError

    def __call__(self, data: pd.DataFrame) -> pd.Series:
        return self._execute(data)

    def __repr__(self):
        # return f"{self.hue} < {self.threshold}"
        raise NotImplementedError


class ConditionGt(ConditionPlot):
    def __init__(self, hue: str, threshold: float):
        super().__init__(hue, threshold)

    def _execute(self, data: pd.DataFrame) -> pd.Series:
        return data[self.hue] > self.threshold

    def __repr__(self):
        return f"{self.hue} > {self.threshold}"


class ConditionLt(ConditionPlot):
    def __init__(self, hue: str, threshold: float):
        super().__init__(hue, threshold)

    def _execute(self, data: pd.DataFrame) -> pd.Series:
        return data[self.hue] < self.threshold

    def __repr__(self):
        return f"{self.hue} < {self.threshold}"


def welch_ttest(data1, data2):
    """Welchのt検定を行う

    Returns:
        _type_: _description_
    """


def plot_distplots(
    data: pd.DataFrame,
    x: str,
    thresholds: List[ConditionPlot],
    xlim: Tuple = (-0.1, 1.1),
    ylim: Tuple = (0.0, 4.0),
):
    aspectes = 5
    fig, axes = plt.subplots(
        len(thresholds), 1, figsize=(aspectes * 4, len(thresholds) * 4)
    )

    for i, threshold in enumerate(thresholds):
        sample = data[threshold(data)]
        sns.kdeplot(sample, x=x, ax=axes[i])
        t, p = stats.ttest_ind(sample[x], data[x], equal_var=False)
        axes[i].set_title(str(threshold) + f" (n={sample.shape[0]}, p={p:.2e})")
        axes[i].set_xlim(xlim)
        axes[i].set_ylim(ylim)  # type: ignore

    fig.tight_layout()
    return fig
