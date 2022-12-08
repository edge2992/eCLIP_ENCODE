# %%
import os
from dotenv import load_dotenv

from src.plot.interaction_metrics.representative import (
    target_report,
    metrics,
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
data = metrics(report)

# %%
data.head()
data.shape

data.sort_values("simpson", ascending=False).reset_index(drop=True).to_csv(
    os.path.join(SAVEDIR, "simpson_sorted_{}.csv".format(BIOSAMPLE)), index=False
)

# %%
