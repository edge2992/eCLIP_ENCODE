# %%
import os

import matplotlib.pyplot as plt
from dotenv import load_dotenv
import numpy as np
import pandas as pd
from typing import Dict

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.similarity_strategy import INTERACTION_SIMILARITY_STRATEGIES
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG
from thesis.utils.compareset import COMPARESET

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(THESIS_FIG_PATH, "results", "protein_interaction")
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results", "protein_interaction")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%

report_data = {
    biosample: {
        k: SampleSetECLIP(condition) for k, condition in COMPARESET[biosample].items()
    }
    for biosample in COMPARESET
}

# %%
d = {}  # type: Dict[str, Dict[str, int]]
for biosample in report_data:
    d[biosample] = {}
    for k, v in report_data[biosample].items():
        print(k)
        print(v.report.shape, v.report["Target label"].nunique())
        d[biosample][k] = v.report["Target label"].nunique()

protein_metrics_table = pd.DataFrame(d).T

protein_metrics_table

# %%

protein_metrics_table.style.format_index(escape="latex", axis=1).to_latex(
    os.path.join(SAVEDIR, "compareset_protein_number.tex"),
    position="tbp",
    position_float="centering",
    hrules=True,
    caption="相互作用するRNA (遺伝子) の種類数で分割した場合のタンパク質の種類数",
    label="tb:compareset_protein",
)
# %% タンパク質の数, 相互作用数、RNAの種類数について
# Gene Jaccard, Gene Simpsonについて

interaction_data = {
    biosample: {
        k: Metrics(SampleSetECLIP(condition).report)(INTERACTION_SIMILARITY_STRATEGIES)
        for k, condition in COMPARESET[biosample].items()
    }
    for biosample in COMPARESET
}

# %%
plot_data: Dict[str, pd.DataFrame] = {
    f"{k} ({biosample})": v
    for biosample in interaction_data
    for k, v in interaction_data[biosample].items()
}  # type: ignore

# %%
plot_data.keys()

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
            os.path.join(SAVEDIR, f"histgram_{'_'.join(metrics_col.split())}_set.pdf"),
            bbox_inches="tight",
            dpi=300,
        )
    return fig


# %%
plot_hist_interaction(
    plot_data, metrics_col="Gene Jaccard", plot_kws=dict(bins=50), save=True
)

# %%
plot_hist_interaction(
    plot_data, metrics_col="Gene Simpson", plot_kws=dict(bins=50), save=True
)

# %% basic data

# %%
meta: Dict[str, pd.DataFrame] = {
    f"{biosample}_{k}": v
    for biosample in interaction_data
    for k, v in interaction_data[biosample].items()
}  # type: ignore

concat_data: pd.DataFrame = pd.concat([value for _, value in meta.items()], keys=meta.keys()).reset_index()  # type: ignore

metrics_table = (
    pd.concat(
        [
            concat_data["level_0"]
            .str.split("_", 1, expand=True)  # type: ignore
            .rename(columns={0: "Biosample", 1: "Condition"}),
            concat_data.loc[:, ["Gene Jaccard", "Gene Simpson"]],
        ],
        axis=1,
    )
    .groupby(["Biosample", "Condition"])
    .agg({"Gene Jaccard": ["count", "mean", "std"], "Gene Simpson": ["mean", "std"]})
)

# %%
metrics_table

# %%

metrics_table.style.format(
    {
        ("Gene Jaccard", "mean"): "{:.2e}",
        ("Gene Jaccard", "std"): "{:.2e}",
        ("Gene Simpson", "mean"): "{:.2e}",
        ("Gene Simpson", "std"): "{:.2e}",
        ("Gene Jaccard", "count"): "{:,d}",
    }
).format_index(escape="latex", axis=1).format_index(escape="latex", axis=0).to_latex(
    os.path.join(SAVEDIR, "compareset_metrics_table.tex"),
    position="tbp",
    position_float="centering",
    hrules=True,
    caption="比較セットごとのタンパク質の相互作用の統計量",
    label="tab:compareset_metrics_basic",
)
# %%
