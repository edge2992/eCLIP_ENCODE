# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.eclip.encodecondition import ECLIP_SAMPLESETS
from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.similarity_strategy import (
    TAPE,
    BlastP,
    DirectStringScore,
    FoldSeekTMScore,
    Jaccard,
    KeywordCosine,
    PeakStrategy,
)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
SAVEDIR = os.path.join(PROJECT_PATH, "src/plot/img", "distplot_metrics")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# %%

STRATEGIES = [
    Jaccard(),
    PeakStrategy(metrics="jaccard"),
    PeakStrategy(metrics="union-intersection"),
    PeakStrategy(metrics="intersection"),
    PeakStrategy(metrics="n_intersections"),
    TAPE(symmetric=True),
    BlastP(symmetric=True),
    KeywordCosine(symmetric=True),
    DirectStringScore(metrics="score", symmetric=True),
    DirectStringScore(metrics="ascore", symmetric=True),
    DirectStringScore(metrics="escore", symmetric=True),
    DirectStringScore(metrics="tscore", symmetric=True),
    FoldSeekTMScore(symmetric=True),
]

# %%

data_list = []
for conditions in ECLIP_SAMPLESETS:
    sampleset = SampleSetECLIP(conditions)
    data = Metrics(sampleset.report)(STRATEGIES)
    data["label"] = str(sampleset)
    data_list.append(data)

# %%

dd = pd.concat(data_list)
target_cols = dd.columns[6:-1]
fig, axes = plt.subplots(len(target_cols), 1, figsize=(5, 5 * len(target_cols)))
for col, ax in zip(target_cols, axes):
    sns.violinplot(data=dd, x=col, y="label", ax=ax)
    ax.set_title(col)
fig.savefig(os.path.join(SAVEDIR, "violinplot_metrics.png"), bbox_inches="tight")


# %%
dd = pd.concat(data_list)
dd.replace(
    {
        "Biosample_name__HepG2_and_Gene_Number_over_1000": "HepG2\nover1000",
        "Biosample_name__HepG2_and_Gene_Number_under_1000": "HepG2\nunder1000",
        "Biosample_name__K562_and_Gene_Number_over_1000": "K562\nover1000",
        "Biosample_name__K562_and_Gene_Number_under_1000": "K562\nunder1000",
    },
    inplace=True,
)
print(dd["Jaccard"].quantile(0.75))
dd["Jaccard_STRONG"] = dd["Jaccard"] > dd["Jaccard"].quantile(0.75)

# %%
target_cols = dd.columns[6:-2]
fig, axes = plt.subplots(len(target_cols), 1, figsize=(5, 5 * len(target_cols)))
for col, ax in zip(target_cols, axes):
    sns.violinplot(data=dd, x="label", y=col, hue="Jaccard_STRONG", ax=ax)
fig.savefig(
    os.path.join(SAVEDIR, "violinplot_metrics_by_jaccard.png"), bbox_inches="tight"
)

# %%
import scipy.stats as stats

target_cols = dd.columns[6:-2]

# for name, grouped in dd.groupby(["label", "Jaccard_STRONG"]):
#     print(name)
#     for col in target_cols:
#         print(col)
#         print(stats.kstest(grouped[col], "norm"))

# %%
results = []
for name, grouped in dd.groupby("label"):
    for col in target_cols:
        for alt in ["two-sided", "less", "greater"]:
            statistic, pvalue = stats.mannwhitneyu(
                grouped[grouped["Jaccard_STRONG"]][col],
                grouped[~grouped["Jaccard_STRONG"]][col],
                alternative=alt,
            )
            results.append(
                {
                    "label": name,
                    "metrics": col,
                    "statistic": statistic,
                    "pvalue": pvalue,
                    "alternative": alt,
                }
            )
# %%
results_df = pd.DataFrame(results).pivot(
    index=["metrics", "label"], columns="alternative", values="pvalue"
)
results_df.head()

# %%
replaced_result = results_df.reset_index().replace(
    {
        "HepG2\nover1000": "HepG2_over1000",
        "HepG2\nunder1000": "HepG2_under1000",
        "K562\nover1000": "K562_over1000",
        "K562\nunder1000": "K562_under1000",
    }
)
for name, grouped in replaced_result.groupby("metrics"):
    fig, ax = plt.subplots()
    sns.heatmap(
        grouped.set_index("label")[["greater", "less", "two-sided"]],
        annot=True,
        ax=ax,
        fmt=".2e",
    )
    ax.set_title(str(name))
    fig.savefig(os.path.join(SAVEDIR, f"heatmap_{name}.png"), bbox_inches="tight")

# %%
