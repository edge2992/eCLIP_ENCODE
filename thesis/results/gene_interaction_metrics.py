# %%
import os

import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import ConditionEq
from src.util.similarity_strategy import Gene_N_Min, Gene_N_Union, Jaccard, Simpson
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
]

interaction_data = {
    k: Metrics(SampleSetECLIP(ConditionEq("Biosample name", k)).report)(
        METRICS, add_description=False  # type: ignore
    )
    for k in ["HepG2", "K562"]
}

# %%
for biosample in interaction_data.keys():
    fig, axes = plt.subplots(1, 2)
    for hue, ax in zip(["Gene Jaccard", "Gene Simpson"], axes):
        sns.scatterplot(
            interaction_data[biosample],
            y="Gene N Union",
            x="Gene N Min",
            ax=ax,
            hue=hue,
            edgecolor="none",
        )
        handlers, labels = ax.get_legend_handles_labels()
        # ax.legend(handlers, labels, loc="upper left", bbox_to_anchor=(1, 1))
        ax.set_xscale("log")
        ax.grid(which="major", ls="-")
        ax.set_xlim(30, 5200)
        ax.set_ylim(-1, 8000)
        # TODO scale
    fig.subplots_adjust(wspace=0.30)
    fig.savefig(
        os.path.join(SAVEDIR, f"interaction_union_min_{biosample}.pdf"),
        bbox_inches="tight",
        dpi=300,
    )

# %%
for biosample in interaction_data.keys():
    print(interaction_data[biosample].max())

# %%
