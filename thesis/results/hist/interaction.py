# %%
import os

import matplotlib.pyplot as plt
from dotenv import load_dotenv
import numpy as np
import pandas as pd
from typing import Dict

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import ConditionEq
from src.util.similarity_strategy import INTERACTION_SIMILARITY_STRATEGIES
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


interaction_data: Dict[str, pd.DataFrame] = {
    k: Metrics(SampleSetECLIP(ConditionEq("Biosample name", k)).report)(
        INTERACTION_SIMILARITY_STRATEGIES, add_description=False
    )
    for k in ["HepG2", "K562"]
}  # type: ignore

# %%
interaction_data["HepG2"].head()

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
        plot_kws = dict(
            bins=np.logspace(np.log10(min_), np.log10(max_), plot_kws["bins"])
        )

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
            os.path.join(SAVEDIR, f"histgram_{'_'.join(metrics_col.split())}.pdf"),
            bbox_inches="tight",
            dpi=300,
        )
    return fig


# %%
for col in interaction_data["HepG2"].columns:
    plot_hist_interaction(
        interaction_data,
        metrics_col=col,
        plot_kws=dict(bins=90),
        x_log=col in ["Peak Intersection", "Peak Jaccard", "Peak N_intersections"],
        save=True,
    )

# %%
