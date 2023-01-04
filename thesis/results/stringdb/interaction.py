# %%
import os
import pandas as pd

import matplotlib.pyplot as plt
from dotenv import load_dotenv
import seaborn as sns
from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    INTERACTION_SIMILARITY_STRATEGIES,
    DirectStringScore,
)
from thesis.utils.reportset import COMPARE_REPORT_SET
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG
from thesis.utils.tex_escape import tex_escape

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results", "stringdb", "interaction")

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %%

for biosample in COMPARE_REPORT_SET:
    for condition in COMPARE_REPORT_SET[biosample]:
        sampleset = COMPARE_REPORT_SET[biosample][condition]

        data: pd.DataFrame = Metrics(sampleset.report)(
            [
                *INTERACTION_SIMILARITY_STRATEGIES,
                DirectStringScore(metrics="score", qcut=True),
            ],
            add_description=False,
        ).replace(
            {
                "STRING Score": {
                    "h": "high-confidence",
                    "m": "midium-confidence",
                    "l": "low-confidence",
                }
            },
        )  # type: ignore

        metrics_table = data.groupby("STRING Score").agg(
            {"Gene Jaccard": ["count", "mean", "std"], "Gene Simpson": ["mean", "std"]}
        )

        metrics_table.style.format(
            {
                ("Gene Jaccard", "mean"): "{:.2e}",
                ("Gene Jaccard", "std"): "{:.2e}",
                ("Gene Simpson", "mean"): "{:.2e}",
                ("Gene Simpson", "std"): "{:.2e}",
                ("Gene Jaccard", "count"): "{:,d}",
            }
        ).format_index(escape="latex", axis=1).format_index(
            escape="latex", axis=0
        ).to_latex(
            os.path.join(
                TB_SAVEDIR,
                f"metrics_table_string_interaction_{biosample}_{condition}.tex",
            ),
            position="htbp",
            position_float="centering",
            hrules=True,
            caption=tex_escape(
                f"STRING Scoreで分けた時の相互作用類似指標の代表値 ({condition} ({biosample}))"
            ),
            label="tab:metrics_string_interaction_{biosample}_{condition}",
        )


# %%
# heavy no meaning
condition = "Protein_GTE_1000"

for biosample in ["HepG2", "K562"]:
    sampleset = COMPARE_REPORT_SET[biosample][condition]
    data: pd.DataFrame = Metrics(sampleset.report)(
        [
            *INTERACTION_SIMILARITY_STRATEGIES,
            DirectStringScore(metrics="score", qcut=True),
        ],
        add_description=False,
    ).replace(
        {
            "STRING Score": {
                "h": "high-confidence",
                "m": "midium-confidence",
                "l": "low-confidence",
            }
        },
    )  # type: ignore

    splot = sns.pairplot(
        data,
        hue="STRING Score",
        vars=["Gene Jaccard", "Gene Simpson", "Peak Jaccard"],
        diag_kind="kde",
        corner=True,
        height=4,
        aspect=1,
        plot_kws=dict(alpha=0.15, edgecolor="none", s=70),
    )

    # splot.fig.savefig(
    #     os.path.join(SAVEDIR, f"interaction_stringdb_metrics_pairplot_{biosample}.pdf"),
    #     dpi=300,
    # )

# %%
