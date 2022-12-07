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
    get_keyword,
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
data = metrics(report)
data.shape

# %%
keyword_experiment_pair: Dict[str, List] = {}

for index, row in tqdm(data.iterrows()):
    keyword1 = get_keyword(row["Dataset_1"])
    keyword2 = get_keyword(row["Dataset_2"])
    intersection_keyword = list(set(keyword1) & set(keyword2))
    for k in intersection_keyword:
        if k not in keyword_experiment_pair:
            keyword_experiment_pair[k] = []
        keyword_experiment_pair[k].append(index)

# %%
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
