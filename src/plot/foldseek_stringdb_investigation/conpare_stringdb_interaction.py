# stringdbのスコアとeCLIPの相互作用のスコアを比較する
# %% import
import os
import textwrap

import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.plot.foldseek_stringdb_investigation.plot_utils import (
    ConditionGt,
    plot_distplots,
)
from src.util.metrics import Metrics
from src.plot.interaction_metrics.representative import target_report
from src.util.similarity_strategy import (
    DirectStringScore,
    Jaccard,
    Simpson,
)

sns.set(font_scale=1.4)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(
    PROJECT_PATH,
    "src/plot/img",
    "foldseek_stringdb_investigation",
    "compare_stringdb_interaction",
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


# %% prepare_plotting data

THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
data: pd.DataFrame = Metrics(report)(
    [
        DirectStringScore(report, metrics="score"),
        DirectStringScore(report, metrics="ascore"),
        DirectStringScore(report, metrics="escore"),
        DirectStringScore(report, metrics="tscore"),
        Simpson(),
        Jaccard(),
    ]
)  # type: ignore


# %%
# keyword_experiment_pair = convert_to_dict_exp_pair_by_keyword(data)
# %% pairplot
data.head()

stringdb_score_nonzero_data = data[data["stringdb_score"] != 0]
print(data.shape, "->", stringdb_score_nonzero_data.shape)

splot = sns.pairplot(
    stringdb_score_nonzero_data,
    x_vars=[
        "stringdb_score",
        "stringdb_ascore",
        "stringdb_escore",
        "stringdb_tscore",
    ],
    y_vars=[
        "Simpson",
        "Jaccard",
    ],
    plot_kws={"alpha": 0.5},
)

description = """stringdb_score vs interaction_score. Among the data obtained in experiments with the HepG2 strain, those with more than 1000 types of interacting genes were used. I also removed the ones where stringdb_score was null."""  # NOQA
splot.fig.suptitle(textwrap.fill(description, width=60), y=1.05, fontsize=16)
splot.fig.subplots_adjust(top=0.8, bottom=0.1, left=0.1, right=0.9)

splot.fig.savefig(
    os.path.join(SAVEDIR, "stringdb_score_interaction_pairplot.png"),
    bbox_inches="tight",
)


# %% distplot
# stringdb_scoreの閾値以上のデータを取り出した時にjaccard, simpsonスコアの分布を確認する 逆もやる

for score_metrics in [
    "stringdb_score",
    "stringdb_ascore",
    "stringdb_escore",
    "stringdb_tscore",
]:
    nonzero_data = data[data[score_metrics] != 0]
    print(data.shape, "->", nonzero_data.shape)
    fig = plot_distplots(
        nonzero_data,
        x="Jaccard",
        thresholds=[ConditionGt(score_metrics, t) for t in [0.1, 0.3, 0.5, 0.7, 0.9]],
        xlim=(-0.1, 0.7),
        ylim=(0.0, 7.0),
    )
    fig.savefig(
        os.path.join(SAVEDIR, f"jaccard_{score_metrics}_distplot.png"),
        bbox_inches="tight",
    )

    fig = plot_distplots(
        nonzero_data,
        x="Simpson",
        thresholds=[ConditionGt(score_metrics, t) for t in [0.1, 0.3, 0.5, 0.7, 0.9]],
        xlim=(-0.1, 0.9),
        ylim=(0.0, 5.0),
    )
    fig.savefig(
        os.path.join(SAVEDIR, f"simpson_{score_metrics}_distplot.png"),
        bbox_inches="tight",
    )
# %%
