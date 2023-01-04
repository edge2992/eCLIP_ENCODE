# proteinの統計量についてGTE1000とGL1000に分けた時の分布を可視化しておく
# %%
import os
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dotenv import load_dotenv

from src.util.similarity_strategy import PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(
    THESIS_FIG_PATH, "results", "protein_interaction", "hist_protein"
)
TB_SAVEDIR = os.path.join(
    THESIS_TB_PATH, "results", "protein_interaction", "hist_protein"
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value


# %%
from thesis.results.protein_interaction.load import data, report_data

d = {}  # type: ignore

for biosample in report_data:
    d[biosample] = {}
    for key in report_data[biosample]:
        d[key] = report_data[biosample][key].proteins

# %%


def plot_hist_interaction(
    data: Dict[str, pd.DataFrame],
    metrics_col: str = "Gene Jaccard",
    plot_kws: Dict = dict(bins=90),
    x_log: bool = False,
    save: bool = False,
):
    plot_data = [data[biosample][metrics_col] for biosample in data]

    if x_log:
        max_ = max([d.max() for d in plot_data])
        min_ = min([d[d != 0].min() for d in plot_data])
        plot_kws["bins"] = np.logspace(np.log10(min_), np.log10(max_), plot_kws["bins"])

    fig, ax = plt.subplots(1, 1, figsize=(16, 4))

    ax.hist(
        plot_data,
        label=list(data.keys()),  # type: ignore
        **plot_kws,
    )
    ax.set_xlabel(metrics_col)
    ax.set_ylabel("Frequency")
    ax.legend()
    if x_log:
        ax.set_xscale("log")
    ax.grid(which="major", ls="-")
    if save:
        fig.savefig(
            os.path.join(SAVEDIR, f"histgram_{'_'.join(metrics_col.split())}_set.pdf"),
            bbox_inches="tight",
            dpi=300,
        )
    return fig


plot_data: Dict[str, pd.DataFrame] = {
    f"{k} ({biosample})": v[[str(s) for s in PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES]]
    for biosample in data
    for k, v in data[biosample].items()
}

for s in PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES:
    plot_hist_interaction(
        plot_data,
        str(s),
        plot_kws=dict(bins=50),
        x_log=str(s) in ["BLASTP Bit avg"],
        save=True,
    )

# %%
