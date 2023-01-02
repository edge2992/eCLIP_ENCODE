# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from src.util.metrics import ProteinMetrics
from src.util.similarity_protein import ProteinSimilarity
from src.util.similarity_strategy import DirectStringScore
from src.util.similarity_strategy.stringdb import STRINGDB_SCORE_METRICS
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

STRINGDB_SCORE_STRATEGIES = [
    DirectStringScore(metrics=met, fillna=None) for met in STRINGDB_SCORE_METRICS
]

data = ProteinMetrics()(STRINGDB_SCORE_STRATEGIES, add_description=True)  # type: ignore
# すでにデータベースから取得した時点で0が入っている
# 逆に参照が失敗したタンパク質にはNanが入っている。6464件

# %%
data.fillna(0, inplace=True)
(data == 0).sum()

# %%
stringdb_metrics_filled = pd.concat(
    [(data != 0).sum(), (data != 0).sum() / data.shape[0] * 100], axis=1  # type: ignore
)
stringdb_metrics_filled.columns = ["Non-zero Entry", "Non-zero Entry Ratio"]
stringdb_metrics_filled.iloc[2:].style.format(
    {"Non-zero Entry": "{:,.0f}", "Non-zero Entry Ratio": "{:.2f}\%"}
).to_latex(
    os.path.join(TB_SAVEDIR, "stringdb_metrics_filled.tex"),
    column_format="lrr",
    position="H",
    position_float="centering",
    hrules=True,
    caption="STRINGDB内の各スコアのレコードの数と割合。",
    label="tab:stringdb_metrics_filled",
)

# %%
# タンパク質ごとに見ると

handler = ProteinSimilarity()
handler.setStrategy(DirectStringScore(metrics="score", fillna=None))
matrix = handler.executeStrategy()
print(matrix.isna().sum().sort_values(ascending=False))

# %%
