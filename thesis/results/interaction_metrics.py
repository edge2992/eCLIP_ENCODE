# %%
import os
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
interaction_data: Dict[str, pd.DataFrame] = {
    k: Metrics(SampleSetECLIP(ConditionEq("Biosample name", k)).report)(
        INTERACTION_SIMILARITY_STRATEGIES, add_description=False
    )
    for k in ["HepG2", "K562"]
}  # type: ignore

# %%


def get_lim(index: str):
    max_ = max(
        [interaction_data[biosample][index].max() for biosample in interaction_data]
    )
    min_ = min(
        [interaction_data[biosample][index].min() for biosample in interaction_data]
    )
    return min_, max_


def plot_jaccard_gene_scatter(
    data: Dict[str, pd.DataFrame],
    condition=lambda df: np.log(df["Peak N_intersections"]),
    colorbar_label: str = "log(peak N_intersections)",
):
    xlim = get_lim("Gene Jaccard")
    ylim = get_lim("Peak Jaccard")
    plot_kws = dict(cmap="Spectral")

    fig, axes = plt.subplots(1, 2)

    for biosample, ax in zip(data, axes):
        mappable = ax.scatter(
            interaction_data[biosample]["Gene Jaccard"],
            interaction_data[biosample]["Peak Jaccard"],
            c=condition(interaction_data[biosample]),  # type: ignore
            **plot_kws,
        )
        fig.colorbar(mappable, ax=ax, label=colorbar_label)
        ax.set_xlabel("Gene Jaccard")
        ax.set_ylabel("Peak Jaccard")
        ax.set_yscale("log")
        ax.set_title(biosample)
        ax.set_xlim(xlim)
        ax.set_ylim(1e-6, ylim[1])
        ax.grid(which="major", ls="-")
    fig.subplots_adjust(wspace=0.3)

    return fig


# %%
fig = plot_jaccard_gene_scatter(
    interaction_data,
    condition=lambda df: np.log(df["Peak N_intersections"]),
    colorbar_label="log(peak N_intersections)",
)
fig.savefig(
    os.path.join(SAVEDIR, "jaccard_peak_gene_scatter_peak_N_intersections.pdf"), dpi=300
)

# %%
fig = plot_jaccard_gene_scatter(
    interaction_data,
    lambda df: df["Peak Union-intersection"],
    "Peak Union-intersection",
)
fig.savefig(os.path.join(SAVEDIR, "jaccard_peak_gene_scatter_peak_union.pdf"), dpi=300)

# %%
interaction_data["HepG2"].head()

# %%


def plot_corr_heatmap(
    data: Dict[str, pd.DataFrame],
    biosample: str = "K562",
    corr_method: str = "spearman",
    save: bool = True,
):
    assert corr_method in ["pearson", "spearman"]
    corr = data[biosample].corr(corr_method)  # type: ignore
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
    ax.set_title(f"{corr_method.capitalize()} Correlation matrix ({biosample})")

    if save:
        fig.savefig(
            os.path.join(SAVEDIR, f"corr_{corr_method}_interaction_{biosample}.pdf"),
            dpi=300,
            bbox_inches="tight",
        )

    return fig


# %%
plot_corr_heatmap(interaction_data, "K562", "pearson")
plot_corr_heatmap(interaction_data, "K562", "spearman")
plot_corr_heatmap(interaction_data, "HepG2", "pearson")
plot_corr_heatmap(interaction_data, "HepG2", "spearman")

# %%
