# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from dotenv import load_dotenv

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import ConditionEq
from src.util.similarity_strategy import INTERACTION_SIMILARITY_STRATEGIES
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(THESIS_FIG_PATH, "results")
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%

set_HepG2 = SampleSetECLIP(ConditionEq("Biosample name", "HepG2"))
data_HepG2 = Metrics(set_HepG2.report)(INTERACTION_SIMILARITY_STRATEGIES)
set_K562 = SampleSetECLIP(ConditionEq("Biosample name", "K562"))
data_K562 = Metrics(set_K562.report)(INTERACTION_SIMILARITY_STRATEGIES)

# %%
data_HepG2.head()

# %%


def get_lim(index: str):
    max_ = max(data_HepG2[index].max(), data_K562[index].max())
    min_ = min(data_HepG2[index].min(), data_K562[index].min())
    return min_, max_


xlim = get_lim("Gene Jaccard")
ylim = get_lim("Peak Jaccard")

fig, axes = plt.subplots(1, 2)

axes[0].scatter(data_HepG2["Gene Jaccard"], data_HepG2["Peak Jaccard"])
axes[0].set_xlabel("Gene Jaccard")
axes[0].set_ylabel("Peak Jaccard")
axes[0].set_yscale("log")
axes[0].set_title("HepG2")
axes[0].set_xlim(xlim)
axes[0].set_ylim(1e-6, ylim[1])

axes[1].scatter(data_K562["Gene Jaccard"], data_K562["Peak Jaccard"])
axes[1].set_xlabel("Gene Jaccard")
axes[1].set_ylabel("Peak Jaccard")
axes[1].set_yscale("log")
axes[1].set_title("K562")
axes[1].set_xlim(xlim)
axes[1].set_ylim(1e-6, ylim[1])

fig.subplots_adjust(wspace=0.3)
fig.savefig(os.path.join(SAVEDIR, "jaccard_peak_gene_scatter.pdf"), dpi=300)


# %%
data_K562.head()

# %%


def plot_corr_heatmap(data, method, biosample, save=True):
    corr = data.iloc[:, 6:].corr(method)
    mask = np.zeros_like(corr)
    mask[np.triu_indices_from(mask)] = True

    fig, ax = plt.subplots(1, 1)

    sns.heatmap(
        corr,
        mask=mask,
        annot=True,
        square=True,
        fmt=".2f",
        vmax=1,
        vmin=-1,
        center=0,
        ax=ax,
    )
    ax.set_title(f"{method.capitalize()} Correlation matrix ({biosample})")

    if save:
        fig.savefig(
            os.path.join(SAVEDIR, f"corr_{method}_interaction_{biosample}.pdf"),
            dpi=300,
            bbox_inches="tight",
        )

    return fig


plot_corr_heatmap(data_K562, "pearson", "K562")
plot_corr_heatmap(data_K562, "spearman", "K562")
plot_corr_heatmap(data_HepG2, "pearson", "HepG2")
plot_corr_heatmap(data_HepG2, "spearman", "HepG2")

# # %%

# %%
