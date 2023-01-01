# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from dotenv import load_dotenv

from src.util.metrics import ProteinMetrics
from src.util.similarity_strategy import (
    TAPE,
    BlastP,
    DirectStringScore,
    FoldSeekTMScore,
    KeywordAA,
)
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


PROTEIN_SIMILARITY_STRATEGIES = [
    BlastP(),
    TAPE(),
    FoldSeekTMScore(),
    DirectStringScore(metrics="score"),
    KeywordAA(),
]
data = ProteinMetrics()(PROTEIN_SIMILARITY_STRATEGIES, add_description=True)
data = data.rename({"FoldSeek TM-Score_avg": "FoldSeek TM-Score"}, axis=1)

# %%
data.describe()
# %%


def plot_corr_heatmap(data, method, save=True):
    corr = data.iloc[:, 2:].corr(method)
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
    ax.set_title(f"{method.capitalize()} Correlation matrix (Protein)")

    if save:
        fig.savefig(
            os.path.join(SAVEDIR, f"corr_{method}_protein.pdf"),
            dpi=300,
            bbox_inches="tight",
        )

    return fig


# %%

plot_corr_heatmap(data, "spearman")
plot_corr_heatmap(data, "pearson")


# %%
