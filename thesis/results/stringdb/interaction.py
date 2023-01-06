# %%
import os
import pandas as pd
import scipy.stats as stats

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

        metrics_table.reindex(
            index=["high-confidence", "midium-confidence", "low-confidence"]
        ).style.format(
            {
                ("Gene Jaccard", "mean"): "{:.2e}",
                ("Gene Jaccard", "std"): "{:.2e}",
                ("Gene Simpson", "mean"): "{:.2e}",
                ("Gene Simpson", "std"): "{:.2e}",
                ("Gene Jaccard", "count"): "{:,d}",
            }
        ).format_index(
            escape="latex", axis=1
        ).format_index(
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
                f"STRING Scoreごとの相互作用類似指標の代表値 ({condition} ({biosample}))"
            ),
            label="tab:metrics_string_interaction_{biosample}_{condition}",
        )
# %%
# 正規分布の検定

stats_result = []
for biosample in COMPARE_REPORT_SET:
    for condition in COMPARE_REPORT_SET[biosample]:
        sampleset = COMPARE_REPORT_SET[biosample][condition]

        data: pd.DataFrame = Metrics(sampleset.report)(  # type: ignore
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
        # ManWhitneyU
        for metrics in ["Gene Jaccard", "Gene Simpson"]:
            for target_label in ["high-confidence", "midium-confidence"]:
                A = data[data["STRING Score"] == target_label][metrics]
                B = data[data["STRING Score"] == "low-confidence"][metrics]
                print(len(A), len(B))
                U1, p = stats.mannwhitneyu(A, B, alternative="greater")
                stats_result.append(
                    {
                        "biosample": biosample,
                        "protein": condition,
                        "A": target_label,
                        "B": "low-confidence",
                        "metrics": metrics,
                        "p-value": p,
                        "U1": U1,
                    }
                )
data_stats = pd.DataFrame(stats_result)
# %%
for biosample in ["HepG2", "K562"]:
    data_stats[data_stats["biosample"] == biosample].drop(
        "biosample", axis=1
    ).style.hide(axis="index").format(
        {"p-value": "{:.2e}", "U1": "{:.2e}"}, escape="latex"
    ).to_latex(
        os.path.join(TB_SAVEDIR, f"mannwhitneyu_{biosample}.tex"),
        position="htbp",
        position_float="centering",
        hrules=True,
        caption=tex_escape(f"Mann-Whitneyの片側U検定 ({biosample})"),
        label=f"tab:mannwhitneyu_{biosample}",
    )
# print(stats.shapiro(data["Gene Jaccard"]))

# %%
# heavy no meaning
condition = "Protein_GTE_1000"

for biosample in ["HepG2", "K562"]:
    sampleset = COMPARE_REPORT_SET[biosample][condition]
    data: pd.DataFrame = Metrics(sampleset.report)(  # type: ignore
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
