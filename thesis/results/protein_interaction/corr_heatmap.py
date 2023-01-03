# %%
import os
from typing import Dict, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    INTERACTION_SIMILARITY_STRATEGIES,
    PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES,
)
from thesis.utils.compareset import COMPARESET
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(
    THESIS_FIG_PATH, "results", "protein_interaction", "corr_heatmap"
)
# TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results", "protein_interaction")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# if not os.path.exists(TB_SAVEDIR):
#     os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%

STRATEGIES = PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES + INTERACTION_SIMILARITY_STRATEGIES


data: Dict[str, Dict[str, pd.DataFrame]] = {
    biosample: {
        key: Metrics(SampleSetECLIP(condition).report)(
            STRATEGIES, add_description=False
        )
        for key, condition in COMPARESET[biosample].items()
    }
    for biosample in COMPARESET.keys()
}  # type: ignore

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


for biosample in data:
    for condition in data[biosample]:
        for corr_method in ["spearman", "pearson"]:
            corr = (
                data[biosample][condition]
                .corr(method=corr_method)  # type: ignore
                .loc[
                    [str(s) for s in PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES],
                    [str(s) for s in INTERACTION_SIMILARITY_STRATEGIES],
                ]
            )
            plot_corr_heatmap(
                corr,
                title=f"{corr_method.capitalize()} Correlation matrix ({condition} ({biosample}))",
                save_path=os.path.join(
                    SAVEDIR, f"{biosample}_{condition}_{corr_method}_corr.pdf"
                ),
            )

# %%
