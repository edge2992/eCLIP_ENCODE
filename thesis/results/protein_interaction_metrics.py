# %%
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv
from typing import Union

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import ConditionEq
from src.util.similarity_strategy import (
    INTERACTION_SIMILARITY_STRATEGIES,
    PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES,
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


def plot_corr_heatmap(
    corr: pd.DataFrame,
    mask: Union[np.ndarray, None] = None,
    title: str = "Correlation matrix",
    save_path: Union[str, None] = None,
):
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
    ax.set_title(title)

    if save_path is not None:
        fig.savefig(
            save_path,
            dpi=300,
            bbox_inches="tight",
        )

    return fig


# %%

STRATEGIES = PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES + INTERACTION_SIMILARITY_STRATEGIES


set_HepG2 = SampleSetECLIP(ConditionEq("Biosample name", "HepG2"))
data_HepG2 = Metrics(set_HepG2.report)(STRATEGIES, add_description=False)
set_K562 = SampleSetECLIP(ConditionEq("Biosample name", "K562"))
data_K562 = Metrics(set_K562.report)(STRATEGIES, add_description=False)

assert isinstance(data_HepG2, pd.DataFrame)
assert isinstance(data_HepG2, pd.DataFrame)

data = {
    "HepG2": data_HepG2,
    "K562": data_K562,
}

# %%
for key in data:
    for corr_method in ["spearman", "pearson"]:
        plot_corr_heatmap(
            data[key]
            .corr(method=corr_method)
            .loc[
                [str(s) for s in PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES],
                [str(s) for s in INTERACTION_SIMILARITY_STRATEGIES],
            ],
            title=f"{corr_method.capitalize()} Correlation matrix ({key}, all)",
            save_path=os.path.join(
                SAVEDIR, f"protein_interaction_{corr_method}_corr_{key}.pdf"
            ),
        )

# %%
