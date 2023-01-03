# %%
import os

import matplotlib.pyplot as plt
from dotenv import load_dotenv
import pandas as pd
from typing import Dict
import numpy as np

from src.util.metrics import ProteinMetrics
from src.util.similarity_strategy import PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(THESIS_FIG_PATH, "results", "hist")
# TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results", "hist")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# if not os.path.exists(TB_SAVEDIR):
#     os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%

protein_data: pd.DataFrame = ProteinMetrics()(
    PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES, add_description=False
)  # type: ignore

# %%
print(protein_data.shape)
protein_data.head()

# %%


def plot_hist_protein(
    data: pd.DataFrame,
    metrics_col: str = "TAPE Cosine",
    plot_kws: Dict = dict(bins=90, rwidth=0.85),
    x_log: bool = False,
    save: bool = False,
):
    plot_data = data[metrics_col]

    if x_log:
        max_ = plot_data.max()
        min_ = plot_data[plot_data != 0].min()
        plot_kws["bins"] = np.logspace(np.log10(min_), np.log10(max_), plot_kws["bins"])

    fig, ax = plt.subplots(1, 1, figsize=(16, 4))
    ax.hist(
        plot_data,
        **plot_kws,
    )
    ax.set_xlabel(metrics_col)
    ax.set_ylabel("Frequency")
    if x_log:
        ax.set_xscale("log")
    ax.grid(which="major", ls="-")
    if save:
        fig.savefig(
            os.path.join(SAVEDIR, f"histgram_{'_'.join(metrics_col.split())}.pdf"),
            bbox_inches="tight",
            dpi=300,
        )
    return fig


for col in protein_data.columns:
    plot_hist_protein(
        protein_data,
        metrics_col=col,
        plot_kws=dict(bins=90, rwidth=0.70, color="tab:orange"),
        x_log=col in ["BLASTP Bit avg"],
        save=True,
    )

# %%
