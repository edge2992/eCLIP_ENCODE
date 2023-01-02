# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.metrics.condition import ConditionEq
from src.util.similarity_strategy import Jaccard, PeakStrategy, Simpson
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(THESIS_FIG_PATH, "appendix")
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "appendix")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%

strategies = [
    Jaccard(),
    Simpson(),
    PeakStrategy(metrics="jaccard"),
    PeakStrategy(metrics="union-intersection"),
    PeakStrategy(metrics="intersection"),
    PeakStrategy(metrics="n_intersections"),
]

set_HepG2 = SampleSetECLIP(ConditionEq("Biosample name", "HepG2"))
data_HepG2 = Metrics(set_HepG2.report)(strategies, add_description=False)
set_K562 = SampleSetECLIP(ConditionEq("Biosample name", "K562"))
data_K562 = Metrics(set_K562.report)(strategies, add_description=False)

# %%
data_HepG2.head()
data_HepG2["BioSample"] = "HepG2"
data_K562["BioSample"] = "K562"
interaction_desc: pd.DataFrame = (
    pd.concat([data_HepG2, data_K562])  # type: ignore
    .rename({"Peak N_intersections": "Peak N\_intersections"}, axis=1)
    .groupby("BioSample")
    .describe()
    .T
)

# %%
idx = pd.IndexSlice
interaction_desc.loc[idx[:, ["mean", "std", "min", "max"]], :].style.format(
    "{:.2e}"
).to_latex(
    os.path.join(TB_SAVEDIR, "interaction_metrics.tex"),
    column_format="lrrrrrr",
    position="H",
    position_float="centering",
    hrules=True,
    caption=("相互作用類似度指標の代表値", "相互作用類似度指標の代表値."),
    label="tab:interaction_metrics",
)
# %%
