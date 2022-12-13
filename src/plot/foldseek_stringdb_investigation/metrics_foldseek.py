# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.util.similarity_protein import ProteinSimilarity
from src.util.similarity_strategy import BlastP, FoldSeekTMScore
from src.util.uniprot import load_uniprot_report
from src.util.similarity_strategy import idmapping_mmcif
from foldseek.wrapper import read_aln_tmscore

sns.set(font_scale=1.4)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "foldseek_stringdb_investigation", "metrics_foldseek"
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# %% 統計値を取得する
similarity = ProteinSimilarity()
similarity.setStrategy(FoldSeekTMScore())
tm_data = similarity.executeStrategy()

similarity.setStrategy(BlastP())
blastp_data = similarity.executeStrategy()

assert all(blastp_data.index == tm_data.index)

# %% foldseekが対称行列ではないらしいので、確認する
N = tm_data.shape[0]
tm_array_list = []
for i in range(N):
    for j in range(i + 1, N):
        tm_array_list.append([tm_data.iloc[i, j], tm_data.iloc[j, i]])

tm_array = np.array(tm_array_list)
fig, ax = plt.subplots()
ax.scatter(x=tm_array[:, 0], y=tm_array[:, 1])
fig.savefig(os.path.join(SAVEDIR, "foldseek_symmetry.png"))

# %% balstpが対称行列ではないらしいので、確認する
N = blastp_data.shape[0]
blastp_array_list = []
for i in range(N):
    for j in range(i + 1, N):
        blastp_array_list.append([blastp_data.iloc[i, j], blastp_data.iloc[j, i]])

blastp_array = np.array(blastp_array_list)

fig, ax = plt.subplots()
ax.scatter(x=blastp_array[:, 0], y=blastp_array[:, 1])
ax.set_xlim(0, 200)
ax.set_ylim(0, 200)
fig.savefig(os.path.join(SAVEDIR, "blastp_symmetry.png"))


# %% foldseekとbitscoreを比較する
concat_data = np.concatenate([tm_array, blastp_array], axis=1)

sample = pd.DataFrame(
    concat_data,
    columns=["foldseek_1", "foldseek_2", "blastp_1", "blastp_2"],
)
sample_nonzero = sample[
    (sample["blastp_1"] < 200) & (sample["blastp_2"] < 200) & (sample != 0).all(axis=1)
].copy()

print(sample.shape, " -> ", sample_nonzero.shape)
sample_nonzero.describe()

# %% bitscoreを閾値とした時のfoldseek_1の分布の推移を見る

THRESHOLDS = range(0, 60, 10)

aspects = 5
fig, axes = plt.subplots(len(THRESHOLDS), 1, figsize=(aspects * 4, 4 * len(THRESHOLDS)))

for i, threshold_bit in enumerate(THRESHOLDS):
    sample = sample_nonzero[
        sample_nonzero[["blastp_1", "blastp_2"]].max(axis=1) >= threshold_bit
    ]
    sns.kdeplot(
        sample,
        x="foldseek_1",
        ax=axes[i],
    )
    axes[i].set_title(f"bit score >= {threshold_bit} (n={sample.shape[0]})")
    axes[i].set_ylim(0, 6.0)
    axes[i].set_xlim(-0.1, 1.1)
fig.tight_layout()
fig.savefig(os.path.join(SAVEDIR, "kdeplot_foldseek_blastp_threshold.png"))
# g = sns.FacetGrid(sample_nonzero, row="label", hue="label", aspect=8, height=2)
# g.map(sns.kdeplot, "foldseek_1", fill=True, alpha=1, lw=1.5, bw_method=0.2)

# %% foldseekとbitscoreの相関を見る

splot = sns.pairplot(data=sample_nonzero, plot_kws={"alpha": 0.2})
splot.fig.suptitle("FoldSeek vs BlastP. Sampled ~10000 pairs")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_foldseek_blastp.png"), bbox_inches="tight"
)
# %% min, max, avgなどの処理をして見る

sample_min_max = pd.DataFrame(
    {
        "foldseek_min": sample_nonzero[["foldseek_1", "foldseek_2"]].min(axis=1),
        "foldseek_max": sample_nonzero[["foldseek_1", "foldseek_2"]].max(axis=1),
        "foldseek_avg": sample_nonzero[["foldseek_1", "foldseek_2"]].sum(axis=1) / 2,
        "blastp_min": sample_nonzero[["blastp_1", "blastp_2"]].min(axis=1),
        "blastp_max": sample_nonzero[["blastp_1", "blastp_2"]].max(axis=1),
        "blastp_avg": sample_nonzero[["blastp_1", "blastp_2"]].sum(axis=1) / 2,
    }
)

splot = sns.pairplot(data=sample_min_max, plot_kws={"alpha": 0.2})
splot.fig.suptitle("FoldSeek vs BlastP. Sampled ~10000 pairs")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_foldseek_blastp_min_max_sum.png"),
    bbox_inches="tight",
)
# %%
data = pd.concat(
    [sample_nonzero[["foldseek_1", "foldseek_2"]], sample_min_max["blastp_avg"]], axis=1
)

fig, ax = plt.subplots()
mappable = ax.scatter(
    x=data["foldseek_1"],
    y=data["foldseek_2"],
    c=data["blastp_avg"],
    alpha=0.2,
    cmap="jet",
)
ax.set_xlabel("foldSeek 1")
ax.set_ylabel("foldSeek 2")
ax.set_title("Foldseek1, 2 coded by BlastP average")
fig.colorbar(mappable)
fig.savefig(
    os.path.join(SAVEDIR, "scatter_foldseek_blastp_avg.png"), bbox_inches="tight"
)
# %% 長さに依存しそう -> タンパク質ごとに見てみる
FOLDSEEK_TM_TSV = os.path.join(PROJECT_PATH, "data/afdb/aln_tmscore_{}.tsv".format(1e9))
foldseek_table = read_aln_tmscore(FOLDSEEK_TM_TSV)
converter = {v: k for k, v in idmapping_mmcif().items()}

foldseek_protein = pd.DataFrame(
    {
        "query": foldseek_table["query"].apply(lambda x: converter[x]),
        "target": foldseek_table["target"].apply(lambda x: converter[x]),
        "TMscore": foldseek_table["TMscore"],
    }
)
fig, ax = plt.subplots(figsize=(100, 10))
sns.boxplot(
    foldseek_protein,
    x="query",
    y="TMscore",
)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
fig.savefig(os.path.join(SAVEDIR, "boxplot_foldseek_protein.png"), bbox_inches="tight")

# %%
metrics_query = foldseek_protein.groupby("query").agg(["mean", "std", "count"])
uniprot_report = load_uniprot_report()

data = pd.merge(
    metrics_query, uniprot_report[["From", "Length"]], left_index=True, right_on="From"
)

fig, ax = plt.subplots()
ax.scatter(data["Length"], data["TMscore", "mean"])
ax.set_xlabel("Length")
ax.set_ylabel("TMscore mean")
fig.savefig(
    os.path.join(SAVEDIR, "scatter_foldseek_protein_length.png"), bbox_inches="tight"
)

# %%
