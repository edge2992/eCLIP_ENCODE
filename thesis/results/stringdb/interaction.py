# %%
import os
import pandas as pd

import matplotlib.pyplot as plt
from dotenv import load_dotenv
import seaborn as sns
from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import ConditionEq
from src.util.similarity_strategy import (
    INTERACTION_SIMILARITY_STRATEGIES,
    DirectStringScore,
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


interaction_data = {
    k: Metrics(SampleSetECLIP(ConditionEq("Biosample name", k)).report)(
        INTERACTION_SIMILARITY_STRATEGIES, add_description=False
    )
    for k in ["HepG2", "K562"]
}

# %%

string_data = {
    k: Metrics(SampleSetECLIP(ConditionEq("Biosample name", k)).report)(
        DirectStringScore(metrics="score"), add_description=False
    )
    for k in ["HepG2", "K562"]
}

# %%
idx = pd.IndexSlice
for biosample in ["HepG2", "K562"]:
    desc = (
        pd.concat(
            [
                interaction_data[biosample].rename(
                    {"Peak N_intersections": "Peak N\_intersections"}, axis=1
                ),
                pd.cut(
                    string_data[biosample],  # type: ignore
                    bins=[-0.1, 0.4, 0.6, 1],
                    labels=["low-confidence", "midium-confidence", "high-confidence"],
                ),
            ],
            axis=1,
        )
        .groupby("STRING Score")
        .describe()
        .T
    )
    desc.loc[idx[:, ["count", "mean", "std"]], :].style.format(
        {"mean": "{:.2f}", "std": "{:.2f}", "count": "{:,.0f}"}
    ).to_latex(
        os.path.join(
            TB_SAVEDIR, f"interaction_metrics_grouped_stringdb_{biosample}.tex"
        ),
        position="htbp",
        position_float="centering",
        hrules=True,
        caption=f"相互作用類似度指標の代表値 ({biosample}",
        label=f"tab:interaction_metrics_grouped_stringdb_{biosample}",
    )


# %%
# heavy no meaning
for biosample in ["HepG2", "K562"]:
    data = pd.concat(
        [
            interaction_data[biosample],
            pd.cut(
                string_data[biosample],  # type: ignore
                bins=[-0.1, 0.4, 0.6, 1],
                labels=["low-confidence", "midium-confidence", "high-confidence"],
            ),
        ],
        axis=1,
    )

    splot = sns.pairplot(
        data,
        hue="STRING Score",
        vars=["Gene Jaccard", "Gene Simpson", "Peak Jaccard"],
        diag_kind="kde",
        corner=True,
        height=4,
        aspect=1,
        plot_kws=dict(alpha=0.15, edgecolor="none", s=70),
    )

    splot.fig.savefig(
        os.path.join(SAVEDIR, f"interaction_stringdb_metrics_pairplot_{biosample}.pdf"),
        dpi=300,
    )

# %%
