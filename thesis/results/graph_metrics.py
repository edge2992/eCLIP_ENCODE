# %% import packages
import os

import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from src.eclip.dataset import Dataset
from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics.condition import ConditionEq
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

# %% peak_gene_scatterplot


def make_edges(report: pd.DataFrame) -> pd.DataFrame:
    d = []
    for _, row in report.iterrows():
        dataset = Dataset(row)
        d.append({"Protein": dataset.protein, "RNA": dataset.genes})
    return pd.DataFrame(d).explode("RNA").reset_index()


# %%

set_HepG2 = SampleSetECLIP(ConditionEq("Biosample name", "HepG2"))
set_K562 = SampleSetECLIP(ConditionEq("Biosample name", "K562"))

edges_HepG2 = make_edges(set_HepG2.report)
edges_K562 = make_edges(set_K562.report)

# %%


def get_scale(from_, to_):
    HepG2_degree = edges_HepG2.groupby(from_).nunique()[to_]
    K562_degree = edges_K562.groupby(from_).nunique()[to_]

    max_degree = max(HepG2_degree.max(), K562_degree.max())
    min_degree = min(HepG2_degree.min(), K562_degree.min())

    return list(range(min_degree, max_degree + 1))


def make_degree_distribution(degrees: pd.Series, ax, xs=None):
    n = len(degrees)
    if xs is None:
        max_degree = degrees.max()
        min_degree = degrees.min()
        xs = list(range(min_degree, max_degree + 1))
    ys = [len([degree for degree in degrees if degree >= x]) / n for x in xs]

    ax.scatter(xs, ys, s=60)
    ax.set_xlabel("Degree (log)")
    ax.set_ylabel("Frequency of Degree (log)")
    ax.set_xlim(xs[0], xs[-1] + 1)
    ax.set_yscale("log")
    ax.set_xscale("log")
    ax.grid(which="major", ls="-")


def make_Protein_degree_distribution(edges: pd.DataFrame, ax, xs=None):
    degrees = edges.groupby("Protein").nunique()["RNA"]
    make_degree_distribution(degrees, ax, xs)


def make_RNA_degree_distribution(edges: pd.DataFrame, ax, xs=None):
    degrees = edges.groupby("RNA").nunique()["Protein"]
    make_degree_distribution(degrees, ax, xs)


# %%
fig, ax = plt.subplots(1, 1)
xs = get_scale("RNA", "Protein")
make_RNA_degree_distribution(edges_HepG2, ax, xs)
make_RNA_degree_distribution(edges_K562, ax, xs)
ax.set_title("Degree distribution of RNA.")
ax.legend(["HepG2", "K562"])
fig.savefig(os.path.join(SAVEDIR, "RNA_degree_distribution.pdf"), dpi=300)


# %%

fig, ax = plt.subplots(1, 1)
xs = get_scale("Protein", "RNA")
make_Protein_degree_distribution(edges_HepG2, ax, xs)
make_Protein_degree_distribution(edges_K562, ax, xs)
ax.set_title("Degree distribution of Protein.")
ax.legend(["HepG2", "K562"])
fig.savefig(os.path.join(SAVEDIR, "Protein_degree_distribution.pdf"), dpi=300)


# %%
fig, axes = plt.subplots(1, 2)

xs = get_scale("Protein", "RNA")
make_Protein_degree_distribution(edges_HepG2, axes[0], xs)
make_Protein_degree_distribution(edges_K562, axes[0], xs)
axes[0].set_xlabel("Degree of Protein Node (log)")
axes[0].legend(["HepG2", "K562"], loc="lower left")

xs = get_scale("RNA", "Protein")
make_RNA_degree_distribution(edges_HepG2, axes[1], xs)
make_RNA_degree_distribution(edges_K562, axes[1], xs)
axes[1].set_xlabel("Degree of RNA Node (log)")
axes[1].legend(["HepG2", "K562"], loc="lower left")
fig.subplots_adjust(wspace=0.3)
fig.savefig(os.path.join(SAVEDIR, "degree_distribution.pdf"), dpi=300)


# %%
def graph_metrics(edges: pd.DataFrame):
    return {
        "Edges": len(edges),
        "Protein Nodes": edges["Protein"].nunique(),
        "RNA Nodes": edges["RNA"].nunique(),
    }


# %%
pd.DataFrame(
    {
        "HepG2": graph_metrics(edges_HepG2),
        "K562": graph_metrics(edges_K562),
    }
).to_latex(
    os.path.join(TB_SAVEDIR, "graph_metrics.tex"),
    column_format="lrr",
    caption="Graph metrics for HepG2 and K562.",
    label="tab:graph_metrics",
)

# %%
