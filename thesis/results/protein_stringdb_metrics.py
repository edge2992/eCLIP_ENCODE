# %%
import os

import pandas as pd
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


METRICS = PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES + STRINGDB_SCORE_STRATEGIES
data = ProteinMetrics()(METRICS, add_description=False).fillna(0)

# %%
data.describe()
# %%

# STRING SCORE
# score < 0.4 low-confidence
# 0.4 <= score < 0.6 medium-confidence
# score >= 0.6 high-confidence

# %%
d = {}
for label in [str(s) for s in STRINGDB_SCORE_STRATEGIES]:
    d[label] = pd.cut(
        data[label],
        bins=[-0.1, 0.4, 0.6, 1],
        labels=["low-confidence", "midium-confidence", "high-confidence"],
    )

df_string_label = pd.DataFrame(d)

# %%
pd.concat(
    [
        data.iloc[:, : len(PROTEIN_SIMILARITY_SYMMETRIC_STRATEGIES)],  # type: ignore
        df_string_label,
    ],
    axis=1,
)

# %%
val_count_list = []
for col in df_string_label.columns:
    val_count_list.append(df_string_label[col].value_counts())

pd.DataFrame(val_count_list)[
    ["high-confidence", "midium-confidence", "low-confidence"]
].style.format("{:,d}").to_latex(
    os.path.join(TB_SAVEDIR, "protein_stringdb_label_metrics.tex"),
    column_format="lrrr",
    position="H",
    position_float="centering",
    hrules=True,
    caption="Protein STRINGDB metrics label",
    label="tab:protein_stringdb_metrics_label",
)


# %%
