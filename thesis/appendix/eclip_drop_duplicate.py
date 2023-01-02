# %%
import os

import matplotlib.pyplot as plt
from dotenv import load_dotenv
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

# %%

reports = {
    k: SampleSetECLIP(ConditionEq("Biosample name", k)).report
    for k in ["HepG2", "K562", "adrenal gland"]
}

# %%
for biosample in ["HepG2", "K562", "adrenal gland"]:
    data = reports[biosample]
    print(biosample)
    print(
        data[data.duplicated("Target label")][
            ["Accession", "Target label"]
        ].sort_values("Target label")
    )

# %%
for biosample in ["HepG2", "K562", "adrenal gland"]:
    print(biosample)
    print(
        reports[biosample]
        .sort_values("Date created")["Date created"]
        .map(lambda x: x.split("-")[0])
        .value_counts()
    )

# %%
