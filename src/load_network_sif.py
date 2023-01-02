# eCLIPのDatasetをCytoscapeで読み込めるようにsif形式で出力する


# node_atributes
# %%
import os

import pandas as pd
import networkx as nx
from dotenv import load_dotenv

from src.eclip.dataset import Dataset
from src.eclip.encodecondition import ECLIP_SAMPLESETS
from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.uniprot import load_uniprot_report
from src.util.similarity_strategy import (
    PeakStrategy,
    TAPE,
    BlastP,
    DirectStringScore,
    FoldSeekTMScore,
    Jaccard,
    KeywordCosine,
)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(PROJECT_PATH, "data", "network")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# %%
sampleset = SampleSetECLIP(ECLIP_SAMPLESETS[0])

# %%
NODE_ATTRIBUTES = [
    "ID",
    "Title",
    "Accession",
    "Dataset",
    "Biosample name",
    "Target label",
]


df_nodes = sampleset.report[NODE_ATTRIBUTES]

# %%
keywords = sampleset.report.apply(lambda row: " ".join(Dataset(row).keywords), axis=1)
keywords.name = "Keywords"

gene_num = sampleset.report.apply(lambda row: len(Dataset(row).genes), axis=1)
gene_num.name = "Interacting gene count"

uniprot = load_uniprot_report()

# uniprotの情報を追加
df_nodes = pd.merge(
    pd.concat([sampleset.report[NODE_ATTRIBUTES], keywords, gene_num], axis=1),
    uniprot,
    left_on="Target label",
    right_on="From",
    how="left",
)

# %%

# %%
STRATEGIES = [
    Jaccard(),
    PeakStrategy(metrics="jaccard"),
    PeakStrategy(metrics="union-intersection"),
    PeakStrategy(metrics="intersection"),
    PeakStrategy(metrics="n_intersections"),
    TAPE(symmetric_method="avg"),
    BlastP(symmetric_method="avg"),
    KeywordCosine(),
    DirectStringScore(metrics="score"),
    DirectStringScore(metrics="ascore"),
    DirectStringScore(metrics="escore"),
    DirectStringScore(metrics="tscore"),
    FoldSeekTMScore(symmetric_method="avg"),
]

df_edges = Metrics(sampleset.report)(
    STRATEGIES, add_description=True, report_columns=["ID"]
)

# %%
# %%
assert isinstance(df_edges, pd.DataFrame)
df_edges.rename(columns=lambda x: x.replace(" ", "_").replace("-", "_"), inplace=True)
df_nodes.rename(columns=lambda x: x.replace(" ", "_").replace("-", "_"), inplace=True)

G = nx.from_pandas_edgelist(
    df_edges, source="ID_1", target="ID_2", edge_attr=list(df_edges.columns[2:])
)

node_attr = df_nodes.set_index("ID").to_dict(orient="index")
nx.set_node_attributes(G, node_attr)

# %%
nx.write_gml(
    G,
    os.path.join(SAVEDIR, "eclip_HepG2_over1000.gml"),
)

# %%
