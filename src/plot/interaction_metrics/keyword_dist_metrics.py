# 調査対象の実験ペアについて、キーワードごとに統計情報の分布を可視化する
# TAPE, blastp, keyword, simpsonの分布を可視化する

# %%
import os
from dotenv import load_dotenv
from typing import List, Dict
import pandas as pd
from tqdm import tqdm
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from src.plot.interaction_metrics.representative import (
    similarity_strategy_dict,
    convert_to_dict_exp_pair_by_keyword,
    target_report,
    metrics,
)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
# %%

SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "interaction_metrics", "keyword_dist_metrics"
)
THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"
# BIOSAMPLE = "K562"
TOPN = 200

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


def plot_pairplot_by_keyword(
    data: pd.DataFrame, keyword: str, keyword_exp_pair: Dict[str, List]
):
    # 統計情報をkeywordごとにグラフで表示する
    # TAPE, blastp, keyword, simpsonの分布を可視化する
    def make_hue_label():
        sample = data.copy()
        sample["label"] = "other"
        index_list: List[int] = keyword_exp_pair[keyword]
        sample.iloc[index_list, sample.columns.get_loc("label")] = keyword  # type: ignore
        return sample

    sample = make_hue_label()
    sample["log_blastp (bit score)"] = np.log1p(sample["blastp"])

    palette = {
        keyword: sns.color_palette()[1],  # type: ignore
        "other": sns.color_palette()[0],  # type: ignore
    }

    splot = sns.pairplot(
        sample,
        vars=["TAPE", "keyword", "log_blastp (bit score)", "simpson"],
        hue="label",
        corner=True,
        palette=palette,
        plot_kws={"alpha": 0.3},
    )

    splot.fig.suptitle(
        "pairplot highlighted by [{}, #{}].".format(
            keyword, len(keyword_exp_pair[keyword])
        )
    )
    splot.fig.subplots_adjust(top=0.95)

    return splot


# %%

report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
data = metrics(report, *similarity_strategy_dict())

keyword_experiment_pair = convert_to_dict_exp_pair_by_keyword(data)
keyword_num_series = pd.Series(
    {key: len(values) for key, values in keyword_experiment_pair.items()}
).sort_values(ascending=False)
# %%

for index, (keyword, num) in enumerate(keyword_num_series.items()):
    if num < 10:
        continue
    # 統計情報をグラフで表示する
    # TAPE, blastp, keyword, simpsonの分布を可視化する
    keyword = str(keyword)
    splot = plot_pairplot_by_keyword(data, keyword, keyword_experiment_pair)
    splot.fig.savefig(
        os.path.join(
            SAVEDIR,
            "{}_pairplot_highlighted_by_{}.png".format(
                index, "_".join(keyword.split())
            ),
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.clf()
    plt.close()

# %%


def calc_precision_at_k(
    data: pd.DataFrame, keyword_exp_pair: Dict[str, List[int]], keyword: str, k: int
):
    """simpsonスコアの上位K個を正解とする場合
    keywordを持つものは上位K個のうち何割入っているか"""
    if k == 0 or keyword not in keyword_exp_pair:
        return 0.0
    true_k = data.sort_values("simpson", ascending=False).iloc[:k].index.to_list()
    precision_at_k = len(set(keyword_exp_pair[keyword]) & set(true_k)) / k
    return precision_at_k


def recall_at_k(
    data: pd.DataFrame, keyword_exp_pair: Dict[str, List[int]], keyword: str, k: int
):
    """simpsonスコアの上位K個を正解とする場合
    keywordを持つものはkeywordをもつもの全体のうち何割入っているか"""
    if k == 0 or keyword not in keyword_exp_pair:
        return 0.0
    true_k = data.sort_values("simpson", ascending=False).iloc[:k].index.to_list()
    recall_at_k = len(set(keyword_exp_pair[keyword]) & set(true_k)) / len(
        keyword_exp_pair[keyword]
    )
    return recall_at_k


# %%
result_metrics_list = []
for keyword, num in tqdm(keyword_num_series.items()):
    if num < 10:
        continue
    for k in range(0, 300, 5):
        keyword = str(keyword)
        result_metrics_list.append(
            {
                "keyword": keyword,
                "precision@K": calc_precision_at_k(
                    data, keyword_experiment_pair, keyword, k
                ),
                "recall@K": recall_at_k(data, keyword_experiment_pair, keyword, k),
                "K": k,
            }
        )
result_metrics = pd.DataFrame(result_metrics_list)
result_metrics.head()

# %%
FIG_N = 3
fig, axes = plt.subplots(FIG_N, 1, figsize=(10, FIG_N * 5))
print(axes)

for ax, keys in zip(axes, np.array_split(result_metrics["keyword"].unique(), FIG_N)):
    for keyword in keys:
        sample = result_metrics[result_metrics["keyword"] == keyword]
        ax.plot(sample["recall@K"], sample["precision@K"], label=keyword)
    ax.set_ylabel("Precision@K")
    ax.set_xlabel("Recall@K")
    ax.set_xlim(0, 0.6)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.5), ncol=1)

fig.savefig(
    os.path.join(SAVEDIR, "precision_recall_curve.png"),
    dpi=300,
    bbox_inches="tight",
)
# %%
