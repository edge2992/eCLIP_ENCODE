# %%
import os

import pandas as pd
import matplotlib.pyplot as plt
from dotenv import load_dotenv

from src.util.metrics import ProteinMetrics
from src.util.similarity_strategy import (
    STRINGDB_SCORE_STRATEGIES,
)
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
SAVEDIR = os.path.join(THESIS_FIG_PATH, "proposed")
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "proposed")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%
data = ProteinMetrics()(
    [string_strategy.set_qcut(True) for string_strategy in STRINGDB_SCORE_STRATEGIES],
    add_description=False,
)

# %%
pd.DataFrame([data[col].value_counts() for col in data.columns], dtype=int).fillna(0)[
    ["high-confidence", "midium-confidence", "low-confidence"]
].astype(int).style.format("{:,d}").to_latex(
    os.path.join(TB_SAVEDIR, "protein_stringdb_categorical_metrics.tex"),
    column_format="lrrr",
    position="tbp",
    position_float="centering",
    hrules=True,
    caption="STRINGスコアをカテゴリーに変換した時の各カテゴリーに属するRBPの件数",
    label="tab:stringdb_categorical_metrics",
)

# %%
