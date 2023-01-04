# samplesetごとのSTRINGDBのH, L、Mの数を表形式でまとめる
# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
from dotenv import load_dotenv

from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    INTERACTION_SIMILARITY_STRATEGIES,
    DirectStringScore,
)
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG
from thesis.utils.reportset import COMPARE_REPORT_SET

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]

TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results", "stringdb")

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%
sampledata_num_table = pd.DataFrame(
    index=pd.MultiIndex.from_product(
        [COMPARE_REPORT_SET.keys(), COMPARE_REPORT_SET["HepG2"].keys()]
    ),
    columns=["high-confidence", "midium-confidence", "low-confidence"],
)

for biosample in COMPARE_REPORT_SET:
    for condition in COMPARE_REPORT_SET[biosample]:
        sampleset = COMPARE_REPORT_SET[biosample][condition]
        data = Metrics(sampleset.report)(
            [
                *INTERACTION_SIMILARITY_STRATEGIES,
                DirectStringScore(metrics="score", qcut=True),
            ],
            add_description=False,
        )

        assert isinstance(data, pd.DataFrame)

        sampledata_num_table.loc[(biosample, condition)] = data.replace(
            {
                "STRING Score": {
                    "h": "high-confidence",
                    "m": "midium-confidence",
                    "l": "low-confidence",
                }
            },
            inplace=False,
        )["STRING Score"].value_counts()
# %%

sampledata_num_table.style.format("{:,d}").format_index(
    escape="latex", axis=1
).format_index(escape="latex", axis=0).to_latex(
    os.path.join(TB_SAVEDIR, "stringdb_compareset_count.tex"),
    column_format="llrrr",
    position="htbp",
    position_float="centering",
    hrules=True,
    caption="比較セットごとのタンパク質相互作用のペアの数",
    label="tab:stringdb_compareset_count",
)

# %%
