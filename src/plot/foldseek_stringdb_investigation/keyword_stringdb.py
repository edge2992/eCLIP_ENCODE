# keywordごとにstringdbのスコアの分布を見てみる

# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.plot.interaction_metrics.representative import (
    convert_to_dict_exp_pair_by_keyword,
    metrics,
    target_report,
)
from src.util.similarity_strategy import TAPE, DirectStringScore, Jaccard, Simpson

sns.set(font_scale=1.4)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "foldseek_stringdb_investigation", "keyword_stringdb"
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

    # %%

THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
interaction_strategies = {
    "simpson": lambda report: Simpson(report),
    "jaccard": lambda report: Jaccard(report),
}
protein_strategies = {
    "TAPE": lambda report: TAPE(report),
    "stringdb_score": lambda report: DirectStringScore(report, metrics="score"),
    "stringdb_ascore": lambda report: DirectStringScore(report, metrics="ascore"),
    "stringdb_escore": lambda report: DirectStringScore(report, metrics="escore"),
    "stringdb_tscore": lambda report: DirectStringScore(report, metrics="tscore"),
}

data = metrics(report, protein_strategies, interaction_strategies)

# %%
keyword_experiment_pair = convert_to_dict_exp_pair_by_keyword(data)

# %%


def construct_plotting_data(data: pd.DataFrame, keyword_experiment_pair):
    plot_data = pd.DataFrame()
    for keyword, value in keyword_experiment_pair.items():
        print(keyword, len(value))
        print(data.iloc[value, :].head())
        sample = (
            data.iloc[value, :][
                [
                    "stringdb_score",
                    "stringdb_ascore",
                    "stringdb_escore",
                    "stringdb_tscore",
                ]
            ]
            .copy()
            .reset_index()
        )
        sample["label"] = f"{keyword} (#{len(value)})"
        plot_data = pd.concat([plot_data, sample], axis=0)
    return plot_data


def plot_boxplot_by_keyword(data: pd.DataFrame, metrics: str):
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


# %%

plot_data = construct_plotting_data(data, keyword_experiment_pair)

for met in [
    "stringdb_score",
    "stringdb_ascore",
    "stringdb_escore",
    "stringdb_tscore",
]:
    fig = plot_boxplot_by_keyword(plot_data, met)
    fig.savefig(os.path.join(SAVEDIR, f"boxplot_{met}.png"), bbox_inches="tight")

# %%
