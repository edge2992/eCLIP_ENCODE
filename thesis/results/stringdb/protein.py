# %%
import os

import matplotlib.pyplot as plt
from dotenv import load_dotenv
import seaborn as sns
import numpy as np

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
data["BLASTP Bit avg"] = np.log(data["BLASTP Bit avg"])
# %%
data.describe()
data.head()
# %%

# heavy no meaning

splot = sns.pairplot(
    data.rename(
        {
            "BLASTP Bit avg": "BLASTP Bit (log)",
            "FoldSeek TM-Score avg": "FoldSeek TM-Score",
        },
        axis=1,
    ),
    vars=["BLASTP Bit (log)", "TAPE Cosine", "FoldSeek TM-Score", "Keyword AA"],
    hue="STRING Score",
    diag_kind="kde",
    corner=True,
    height=4,
    aspect=1,
    plot_kws=dict(alpha=0.15, edgecolor="none", s=70),
)

splot.fig.savefig(
    os.path.join(SAVEDIR, "protein_stringdb_metrics_pairplot.pdf"), dpi=300
)

# %%
