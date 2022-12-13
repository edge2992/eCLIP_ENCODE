from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


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
