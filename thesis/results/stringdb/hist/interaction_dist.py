# %%

import os
import pandas as pd
from typing import Dict
import numpy as np

import matplotlib.pyplot as plt
from dotenv import load_dotenv
from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    INTERACTION_SIMILARITY_STRATEGIES,
    DirectStringScore,
)
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG
from thesis.utils.reportset import COMPARE_REPORT_SET

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(THESIS_FIG_PATH, "results", "stringdb")
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results", "stringdb")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%

import seaborn as sns


def plot_hist_interaction_stringdb(
    data: pd.DataFrame,
    metrics_col: str = "Gene Jaccard",
    plot_kws: Dict = dict(bins=90),
    pdf_suffix: str = "",
    x_log: bool = False,
    save: bool = False,
):
    plot_data = {
        key: grouped[metrics_col] for key, grouped in data.groupby("STRING Score")
    }
    color_dict = {
        "high-confidence": "tab:blue",
        "midium-confidence": "tab:orange",
        "low-confidence": "tab:green",
    }

    if x_log:
        max_ = max([d.max() for d in plot_data.values()])
        min_ = min([d[d != 0].min() for d in plot_data.values()])
        plot_kws["bins"] = np.logspace(np.log10(min_), np.log10(max_), plot_kws["bins"])

    fig, ax = plt.subplots(1, 1, figsize=(16, 4))
    ax.hist(
        plot_data.values(),  # type: ignore
        label=list(plot_data.keys()),  # type: ignore
        **plot_kws,
    )
    ax.legend()
    for key, dd in plot_data.items():
        sns.kdeplot(data=dd, ax=ax, label=key, lw=3, color=color_dict[str(key)])
    ax.set_xlabel(metrics_col)
    ax.set_ylabel("Density")
    if x_log:
        ax.set_xscale("log")
    ax.grid(which="major", ls="-")
    if save:
        fig.savefig(
            os.path.join(
                SAVEDIR,
                f"density_histgram_{'_'.join(metrics_col.split())}_{pdf_suffix}.pdf",
            ),
            bbox_inches="tight",
            dpi=300,
        )
    return fig


# %%
for biosample in COMPARE_REPORT_SET:
    for condition in COMPARE_REPORT_SET[biosample]:
        sampleset = COMPARE_REPORT_SET[biosample][condition]
        data = Metrics(sampleset.report)(
            [
                *INTERACTION_SIMILARITY_STRATEGIES,
                DirectStringScore(metrics="score", qcut=True),
            ],
            add_description=False,
        ).replace(
            {
                "STRING Score": {
                    "h": "high-confidence",
                    "m": "midium-confidence",
                    "l": "low-confidence",
                }
            },
        )

        assert isinstance(data, pd.DataFrame)

        for strategies in INTERACTION_SIMILARITY_STRATEGIES:

            fig = plot_hist_interaction_stringdb(
                data,
                metrics_col=str(strategies),
                plot_kws=dict(bins=30, density=True),
                x_log=str(strategies)
                in ["Peak Jaccard", "Peak Intersection", "Peak N_intersections"],
                pdf_suffix=f"{biosample}_{condition}",
                # save=True,
                save=False,
            )
            plt.close()


# %%
