# %%
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.util.similarity_protein import ProteinSimilarity
from src.util.similarity_strategy import FoldSeekTMScore, BlastP

sns.set(font_scale=1.4)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "foldseek_stringdb_investigation", "metrics_foldseek"
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# %%
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
    concat_data[np.random.choice(concat_data.shape[0], 10000, replace=False), :],
    columns=["foldseek_1", "foldseek_2", "blastp_1", "blastp_2"],
)
sample = sample[(sample["blastp_1"] < 200) & (sample["blastp_2"] < 200)]

splot = sns.pairplot(data=sample, plot_kws={"alpha": 0.2})
splot.fig.suptitle("FoldSeek vs BlastP. Sampled ~10000 pairs")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_foldseek_blastp.png"), bbox_inches="tight"
)
# %%

sample = pd.DataFrame(
    {
        "foldseek_max": np.max(tm_array, axis=1),
        "foldseek_min": np.min(tm_array, axis=1),
        "foldseek_mean": np.mean(tm_array, axis=1),
        "blastp_max": np.max(blastp_array, axis=1),
        "blastp_min": np.min(blastp_array, axis=1),
        "blastp_mean": np.mean(blastp_array, axis=1),
    }
)
sample = sample[(sample["blastp_min"] < 200) & (sample["blastp_max"] < 200)]

splot = sns.pairplot(data=sample.sample(10000), plot_kws={"alpha": 0.2})
splot.fig.suptitle("FoldSeek vs BlastP. Sampled ~10000 pairs")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_foldseek_blastp_min_max_sum.png"),
    bbox_inches="tight",
)
# %%
