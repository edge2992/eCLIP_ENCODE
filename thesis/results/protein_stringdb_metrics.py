# %%
import os

import matplotlib.pyplot as plt
from dotenv import load_dotenv

from src.util.metrics import ProteinMetrics
from src.util.similarity_strategy import (
    PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES,
    STRINGDB_SCORE_STRATEGIES,
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


METRICS = PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES + [
    st.set_qcut(True) for st in STRINGDB_SCORE_STRATEGIES
]
data = ProteinMetrics()(METRICS, add_description=False)

# %%
data.describe()
data.head()
# %%

# STRING SCORE
# score < 0.4 low-confidence
# 0.4 <= score < 0.6 medium-confidence
# score >= 0.6 high-confidence
# %%
