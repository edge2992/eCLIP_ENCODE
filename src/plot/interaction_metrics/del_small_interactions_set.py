# 相互作用の少ない実験を排除する
# 優先度高い


# %%
import os
import numpy as np
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv
from scipy.spatial.distance import squareform

from src.plot.util.process_report import count_gene
from src.util.bedfile import load_replicateIDR_report
from src.util.similarity_protein import InteractionSimilarity, ProteinSimilarity
from src.util.similarity_strategy import TAPE, BlastP, KeywordCosine, Simpson

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
# %%

SAVEDIR = os.path.join(PROJECT_PATH, "src/plot/img", "del_small_interactions_set")
THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"
TOPN = 200

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


def target_report(threshold_gene_num: int, biosample: str):
    report = load_replicateIDR_report()
    report = report[report["Biosample name"] == biosample].reset_index(drop=True)
    target_report = report[
        (count_gene(report, lambda row: row["Dataset"]) >= threshold_gene_num).to_list()
    ]
    return target_report


def metrics(report: pd.DataFrame):
    p_similairty = ProteinSimilarity()
    inter_similarity = InteractionSimilarity()

    def protein_similarity(strategy):
        p_similairty.setStrategy(strategy)
        return p_similairty.flatten_tri(
            p_similairty.executeStrategy(transform=True), False
        )

    def interaction_similarity(strategy):
        inter_similarity.setStrategy(strategy)
        return inter_similarity.flatten_tri(inter_similarity.executeStrategy(), False)

    REPORT_COLUMNS = ["Dataset", "Target label", "Biosample name"]

    data = pd.DataFrame(
        {
            "TAPE": protein_similarity(TAPE(report=report)),
            "keyword": protein_similarity(KeywordCosine(report=report)),
            "blastp": protein_similarity(BlastP(report=report)),
            "simpson": interaction_similarity(Simpson(report=report)),
        }
    )
    index_n = np.where(np.triu(squareform(np.ones(data.shape[0]))))
    desc = pd.concat(
        [
            report.iloc[index_n[0]]
            .loc[:, REPORT_COLUMNS]
            .reset_index(drop=True)
            .add_suffix("_1"),
            report.iloc[index_n[1]]
            .loc[:, REPORT_COLUMNS]
            .reset_index(drop=True)
            .add_suffix("_2"),
        ],
        axis=1,
    )
    return pd.concat([desc, data], axis=1).sort_values("TAPE").reset_index(drop=True)


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
data = metrics(report)
# %%
tape_data = data.sort_values("TAPE").reset_index(drop=True)
tape_data["label"] = [
    "TOP{}".format(TOPN) if i else "others" for i in data.index < TOPN
]

splot = sns.pairplot(
    tape_data,
    hue="label",
    vars=["TAPE", "keyword", "blastp", "simpson"],
    corner=True,
    palette={
        "TOP{}".format(TOPN): sns.color_palette()[1],  # type: ignore
        "others": sns.color_palette()[0],  # type: ignore
    },
    plot_kws={"alpha": 0.3},
)
splot.fig.suptitle("pairplot of TAPE, keyword, blastp, simpson sorted by TAPE")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_sorted_by_TAPE.png"), dpi=300, bbox_inches="tight"
)

# %%

simpson_data = data.sort_values("simpson", ascending=False).reset_index(drop=True)
simpson_data["label"] = [
    "TOP{}".format(TOPN) if i else "others" for i in data.index < TOPN
]

splot = sns.pairplot(
    simpson_data,
    hue="label",
    vars=["TAPE", "keyword", "blastp", "simpson"],
    corner=True,
    palette={
        "TOP{}".format(TOPN): sns.color_palette()[1],  # type: ignore
        "others": sns.color_palette()[0],  # type: ignore
    },
    plot_kws={"alpha": 0.3},
)
splot.fig.suptitle("pairplot of TAPE, keyword, blastp, simpson sorted by TAPE")
splot.fig.subplots_adjust(top=0.95)
splot.fig.savefig(
    os.path.join(SAVEDIR, "pairplot_sorted_by_simpson.png"),
    dpi=300,
    bbox_inches="tight",
)

# %%
simpson_data.head()

# %%

# simpsonスコアが高い実験ペアの遺伝子のセットを取得して、比較してみる
