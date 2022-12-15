# %%
import os

import pandas as pd
from dotenv import load_dotenv

from src.plot.interaction_metrics.representative import target_report
from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    TAPE,
    BlastP,
    Cosine,
    KeywordCosine,
    Lift,
    Simpson,
)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
# %%

SAVEDIR = os.path.join(PROJECT_PATH, "src/plot/img", "del_small_interactions_set")
THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"
# BIOSAMPLE = "K562"
TOPN = 200

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
data: pd.DataFrame = Metrics(report)(
    [
        TAPE(),
        KeywordCosine(),
        BlastP(symmetric=True, symmetric_method="avg"),
        Simpson(),
        Lift(),
        Cosine(),
    ]
)  # type: ignore
# %%
data.head()
data.shape

data.sort_values("Simpson", ascending=False).reset_index(drop=True).to_csv(
    os.path.join(SAVEDIR, "simpson_sorted_{}.csv".format(BIOSAMPLE)), index=False
)

# %%
