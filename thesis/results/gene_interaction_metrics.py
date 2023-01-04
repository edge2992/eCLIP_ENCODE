# %%
import os
from typing import Dict

import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import ConditionEq
from src.util.similarity_strategy import (
    Gene_N_Max,
    Gene_N_Min,
    Gene_N_Union,
    Jaccard,
    Simpson,
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

METRICS = [
    Jaccard(),
    Simpson(),
    Gene_N_Union(),
    Gene_N_Min(),
    Gene_N_Max(),
]

interaction_data: Dict[str, pd.DataFrame] = {
    k: Metrics(SampleSetECLIP(ConditionEq("Biosample name", k)).report)(
        METRICS, add_description=False  # type: ignore
    )
    for k in ["HepG2", "K562"]
}

# %%


def plot_interaction_gene_n_scatter(
    data: Dict[str, pd.DataFrame],
    biosample: str = "HepG2",
    x: str = "Gene N Min",
    y: str = "Gene N Union",
    save: bool = True,
):
    plot_data: pd.DataFrame = data[biosample]
    plot_kws = dict(cmap="magma")
    fig, axes = plt.subplots(1, 2)
    for hue, ax in zip(["Gene Jaccard", "Gene Simpson"], axes):
        mappable = ax.scatter(plot_data[x], plot_data[y], c=plot_data[hue], **plot_kws)
        fig.colorbar(mappable, ax=ax, label=hue)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_xscale("log")
        ax.grid(which="major", ls="-")
        # for future: 細胞株を両方表示する場合は、軸を揃える
        # ax.set_xlim(30, 5200)
        # ax.set_ylim(-1, 8000)
    fig.subplots_adjust(wspace=0.30)
    if save:
        fig.savefig(
            os.path.join(
                SAVEDIR,
                f"interaction_{'_'.join(x.split())}_{'_'.join(y.split())}_{biosample}.pdf",
            ),
            bbox_inches="tight",
            dpi=300,
        )


# %%
for biosample in interaction_data.keys():
    plot_interaction_gene_n_scatter(
        interaction_data,
        biosample=biosample,
        x="Gene N Min",
        y="Gene N Union",
        save=True,
    )


# %%
for biosample in interaction_data.keys():
    plot_interaction_gene_n_scatter(
        interaction_data,
        biosample=biosample,
        x="Gene N Max",
        y="Gene N Union",
        save=True,
    )
# %%

for biosample in interaction_data.keys():
    plot_interaction_gene_n_scatter(
        interaction_data,
        biosample=biosample,
        x="Gene N Min",
        y="Gene N Max",
        save=True,
    )

# %%
