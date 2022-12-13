# %%
import os
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.util.similarity_protein import ProteinSimilarity
from src.util.similarity_strategy import DirectStringScore, KeywordCosine

sns.set(font_scale=1.4)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "foldseek_stringdb_investigation", "metrics_stringdb"
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


# %%

similarity = ProteinSimilarity()
similarity.setStrategy(DirectStringScore(metrics="score"))
score_data = similarity.executeStrategy()

# %% 対称行列であることを確認する
N = score_data.shape[0]
score_array_list = []
for i in range(N):
    for j in range(i + 1, N):
        score_array_list.append([score_data.iloc[i, j], score_data.iloc[j, i]])

tm_array = np.array(score_array_list)
fig, ax = plt.subplots()
ax.scatter(x=tm_array[:, 0], y=tm_array[:, 1])
# fig.savefig(os.path.join(SAVEDIR, "stringdb_score_symmetry.png"))

# %%

DIRECT_SCORE_METRICS = [
    "score",
    "nscore",
    "fscore",
    "pscore",
    "ascore",
    "escore",
    "dscore",
    "tscore",
]

metrics_dict = {}
similarity = ProteinSimilarity()

for metrics in DIRECT_SCORE_METRICS:
    similarity.setStrategy(DirectStringScore(metrics=metrics))
    metrics_dict[f"stringdb_{metrics}"] = similarity.flatten_tri(
        similarity.executeStrategy(), include_diagonal=False
    )

similarity.setStrategy(KeywordCosine())
metrics_dict["keyword_cosine"] = similarity.flatten_tri(
    similarity.executeStrategy(), include_diagonal=False
)

# %%
data = pd.DataFrame(metrics_dict)
print(data.shape)
(data == 0).sum()

# %%
# https://towardsdatascience.com/sorry-but-sns-distplot-just-isnt-good-enough-this-is-though-ef2ddbf28078


def plot_distplots(
    data: pd.DataFrame,
    x: str,
    hue: str,
    thresholds: List,
    xlim: Tuple = (-0.1, 1.1),
    ylim: Tuple = (0.0, 4.0),
):
    aspectes = 5
    fig, axes = plt.subplots(
        len(thresholds), 1, figsize=(aspectes * 4, len(thresholds) * 4)
    )

    for i, threshold in enumerate(thresholds):
        sample = data[data[hue] < threshold]

        sns.kdeplot(sample, x=x, ax=axes[i])
        axes[i].set_title(f"{hue} < {threshold} (n={sample.shape[0]})")
        axes[i].set_xlim(xlim)
        axes[i].set_ylim(ylim)  # type: ignore

    fig.tight_layout()
    return fig


# %%
data_nonzero = data[data["stringdb_score"] != 0]

fig, ax = plt.subplots()
sns.scatterplot(data_nonzero, x="stringdb_score", y="keyword_cosine", ax=ax, alpha=0.2)
fig.savefig(
    os.path.join(SAVEDIR, "scatterplot_stringdb_score_keyword_cosine.png"),
    bbox_inches="tight",
)

# %%

fig = plot_distplots(
    data_nonzero,
    x="stringdb_score",
    hue="keyword_cosine",
    thresholds=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0],
)
fig.savefig(
    os.path.join(SAVEDIR, "distplots_stringdb_score_keyword_cosine.png"),
    bbox_inches="tight",
)


# %% score同士の関連を調べる
# 調べたいスコア
# score:   combined score
# ascore: coexpression score
# escore: experimental sore
# tscore: textmining score

OBSERVING_METRICS = ["score", "ascore", "escore", "tscore"]

splot = sns.pairplot(
    data_nonzero[[f"stringdb_{metrics}" for metrics in OBSERVING_METRICS]],
    plot_kws={"alpha": 0.2},
)
splot.fig.suptitle(
    "pairplot of stringdb score.\n"
    "stringdb_score: combined score, stringdb_ascore: coexpression score,\n"
    + "stringdb_escore: experimental sore, stringdb_tscore: textmining score",
)
splot.fig.subplots_adjust(top=0.85)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_stringdb_score.png"), bbox_inches="tight"
)


# %%
