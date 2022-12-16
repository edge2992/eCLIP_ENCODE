# 相互作用の少ない実験を排除する
# 優先度高い

# %%
import pandas as pd
import os
import seaborn as sns
from dotenv import load_dotenv

from src.plot.interaction_metrics.representative import (
    target_report,
    describe_dataset_pair,
)
from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    TAPE,
    KeywordCosine,
    BlastP,
    Simpson,
    Lift,
    Cosine,
)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
# %%

SAVEDIR = os.path.join(PROJECT_PATH, "src/plot/img", "del_small_interactions_set")
THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"
# BIOSAMPLE = "K562"
TOPN = 200

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
data: pd.DataFrame = Metrics(report)(
    [
        TAPE(),
        KeywordCosine(),
        BlastP(symmetric=True, symmetric_method="avg"),
        Simpson(),
        Lift(),
        Cosine(),
    ]
)  # type: ignore
# %%
tape_data = data.sort_values("TAPE").reset_index(drop=True)
tape_data["label"] = [
    "TOP{}".format(TOPN) if i else "others" for i in data.index < TOPN
]

splot = sns.pairplot(
    tape_data,
    hue="label",
    vars=["TAPE", "KeywordCosine", "Blastp", "Simpson"],
    corner=True,
    palette={
        "TOP{}".format(TOPN): sns.color_palette()[1],  # type: ignore
        "others": sns.color_palette()[0],  # type: ignore
    },
    plot_kws={"alpha": 0.3},
)
splot.fig.suptitle("pairplot of TAPE, keyword, blastp, simpson sorted by TAPE")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_sorted_by_TAPE_{}.png".format(BIOSAMPLE)),
    dpi=300,
    bbox_inches="tight",
)

# %%

simpson_data = data.sort_values("Simpson", ascending=False).reset_index(drop=True)
simpson_data["label"] = [
    "TOP{}".format(TOPN) if i else "others" for i in data.index < TOPN
]

splot = sns.pairplot(
    simpson_data,
    hue="label",
    vars=["TAPE", "KeywordCosine", "Blastp", "Simpson"],
    corner=True,
    palette={
        "TOP{}".format(TOPN): sns.color_palette()[1],  # type: ignore
        "others": sns.color_palette()[0],  # type: ignore
    },
    plot_kws={"alpha": 0.3},
)
splot.fig.suptitle("pairplot of TAPE, keyword, blastp, simpson sorted by TAPE")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_sorted_by_simpson_{}.png".format(BIOSAMPLE)),
    dpi=300,
    bbox_inches="tight",
)

# %%
simpson_data.head()

# %%
# simpsonスコアが高い実験ペアの遺伝子のセットを取得して、比較してみる
for index, row in simpson_data.head(10).iterrows():
    describe_dataset_pair(row)

# %%
for index, row in tape_data.head(10).iterrows():
    describe_dataset_pair(row)

# %%
