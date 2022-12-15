# 二つの実験のペアごとに遺伝子の種類数を確認しておく

# %%
import os
from typing import List

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.plot.interaction_metrics.representative import (
    convert_to_dict_exp_pair_by_keyword,
    describe_dataset_pair,
    get_keyword,
    target_report,
)
from src.util.metrics import Metrics
from src.util.similarity_protein import InteractionSimilarity
from src.util.similarity_strategy import (
    TAPE,
    BlastP,
    Cosine,
    KeywordCosine,
    Lift,
    Simpson,
)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
# %%

SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "interaction_metrics", "simpson_gene_num"
)
THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"
# BIOSAMPLE = "K562"
TOPN = 200

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
data: pd.DataFrame = Metrics(report)(
    [
        TAPE(),
        KeywordCosine(),
        BlastP(symmetric=True, symmetric_method="avg"),
        Simpson(),
        Lift(),
        Cosine(),
    ]
)  # type: ignore

# %%
inter_similarity = InteractionSimilarity()
inter_similarity.setStrategy(Simpson(report=report))
inter_similarity.executeStrategy()

# %%
inter_similarity.strategy.accession_genes

extra_metrics = []
for index, row in data.iterrows():
    gene1 = inter_similarity.strategy.accession_genes[row["Dataset_1"]]
    gene2 = inter_similarity.strategy.accession_genes[row["Dataset_2"]]
    extra_metrics.append(
        {
            "Gene_1": len(gene1),
            "Gene_2": len(gene2),
            "Gene_1&2": len(set(gene1) & set(gene2)),
        }
    )

extra_data = pd.concat([data, pd.DataFrame(extra_metrics)], axis=1)

# %%
extra_data.head()

# %%

# %%

METRICS_COLUMNS = ["Simpson", "Lift", "Cosine"]
for metric in METRICS_COLUMNS:
    x1 = extra_data["Gene_1"].to_numpy()
    x2 = extra_data["Gene_2"].to_numpy()
    z = extra_data[metric].to_numpy()

    for i in range(len(x1)):
        if x1[i] > x2[i]:
            tmp = x1[i]
            x1[i] = x2[i]
            x2[i] = tmp

    fig, ax = plt.subplots()
    sns.scatterplot(x=x1, y=x2, hue=z, palette="viridis", ax=ax)
    ax.set_xlabel("Gene_1")
    ax.set_ylabel("Gene_2")
    ax.legend(title=metric)

    fig.savefig(os.path.join(SAVEDIR, f"{BIOSAMPLE}_{metric}_gene_num.png"))


# %%

CHECK_TOP_N = 30
keywords: List[str] = []


METRICS_COLUMNS = ["Simpson", "Lift", "Cosine"]
keyword_common = []
for metric in METRICS_COLUMNS:
    for index, row in (
        extra_data.sort_values(metric, ascending=False).head(CHECK_TOP_N).iterrows()
    ):
        describe_dataset_pair(row)
        keyword1 = get_keyword(row["Dataset_1"])
        keyword2 = get_keyword(row["Dataset_2"])
        inter_keyword = set(keyword1) & set(keyword2)
        for k in inter_keyword:
            keyword_common.append({"keyword": k, "metric": metric})

# %%
sample = pd.DataFrame(keyword_common)

vv = sample["keyword"].value_counts()
sample_ex = sample[sample["keyword"].isin(vv[vv > 5].index)]

fig, ax = plt.subplots(figsize=(10, 20))
sns.countplot(sample_ex, y="keyword", hue="metric", ax=ax)
fig.savefig(
    os.path.join(SAVEDIR, f"{BIOSAMPLE}_keyword_common.png"), bbox_inches="tight"
)

# %%
keyword_experiment_pair = convert_to_dict_exp_pair_by_keyword(data)

# %%
splice_data: pd.DataFrame = data.iloc[keyword_experiment_pair["Spliceosome"]]  # type: ignore
splice_data[splice_data["Simpson"] >= 0.75].sort_values(
    "Simpson", ascending=False
).to_csv(os.path.join(SAVEDIR, f"{BIOSAMPLE}_Spliceosome_high_simpson.csv"))


# %%
data[
    (data["TAPE"].rank() > (len(data) / 2))
    & (data["Simpson"].rank(ascending=False) < 300)
].to_csv(os.path.join(SAVEDIR, f"{BIOSAMPLE}_low_TAPE_high_simpson.csv"))

# %%
