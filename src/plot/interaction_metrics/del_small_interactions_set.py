# 相互作用の少ない実験を排除する
# 優先度高い


# %%
import pandas as pd
from src.plot.util.process_report import count_gene
from src.util.bedfile import load_replicateIDR_report

from src.util.similarity_protein import ProteinSimilarity, InteractionSimilarity
from src.util.similarity_strategy import BlastP, TAPE, KeywordCosine, Simpson

# %%
# %%

# %%
similarity = ProteinSimilarity()

similarity.setStrategy(BlastP())
blastp = similarity.executeStrategy(transform=True)

similarity.setStrategy(TAPE())
tape = similarity.executeStrategy(transform=True)

similarity.setStrategy(KeywordCosine())
keyword = similarity.executeStrategy(transform=True)

# %%
similarity = InteractionSimilarity()
similarity.setStrategy(Simpson())
simpson = similarity.executeStrategy()

# %%
simpson.head()

# %%
assert all(tape.columns == simpson.columns)
assert all(tape.columns == keyword.columns)
assert all(tape.columns == blastp.columns)

# %%
THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"


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
data.head()

# %%
data.sort_values("TAPE").head(10)

# %%
index_n = np.where(np.triu(squareform(np.ones(1081))))
desc1 = (
    report.iloc[index_n[0]]
    .loc[:, ["Dataset", "Target label", "Biosample name"]]
    .reset_index(drop=True)
    .add_suffix("_1")
)
desc2 = (
    report.iloc[index_n[1]]
    .loc[:, ["Dataset", "Target label", "Biosample name"]]
    .reset_index(drop=True)
    .add_suffix("_2")
)

# %%
# print(desc2)
# print(desc1)
dd = pd.concat([desc1, desc2], axis=1)
pd.concat([dd, data], axis=1)
# %%
report["Dataset"]


# %%
print(data.head())
data.shape


# %%
data = pd.DataFrame(
    {
        "TAPE": similarity.flatten_tri(
            tape.loc[target_accessions, target_accessions], False
        ),
        "Keyword": similarity.flatten_tri(
            keyword.loc[target_accessions, target_accessions], False
        ),
        "blastp": similarity.flatten_tri(
            blastp.loc[target_accessions, target_accessions], False
        ),
        "simpson": similarity.flatten_tri(
            simpson.loc[target_accessions, target_accessions], False
        ),
    }
)

# %%
from scipy.spatial.distance import squareform
import numpy as np

TOPN = 1000
data["rank"] = data["TAPE"].rank()
mat = np.triu(squareform(data["rank"] <= TOPN))
index_topN = np.where(mat == True)

# %%
protein_data = pd.DataFrame(
    {
        "dataset1": [target_accessions[i] for i in index_topN[0]],
        "dataset2": [target_accessions[i] for i in index_topN[1]],
        "TAPE": squareform(data["TAPE"])[index_topN],
        "blastp": squareform(data["blastp"])[index_topN],
        "keyword": squareform(data["Keyword"])[index_topN],
        "simpson": squareform(data["simpson"])[index_topN],
    }
)
protein_data.head()
# %%
# tapeの順位100位を取って
# dataset1, protien1, cell1, dataset2, protein2, cell2, tape, blastp, keyword, simpsonのデータフレームを作成する

# %%
report = load_replicateIDR_report()[["Dataset", "Target label", "Biosample name"]]
report.head()

# %%
targets = (
    pd.merge(
        pd.merge(
            protein_data, report, left_on="dataset1", right_on="Dataset", how="inner"
        ),
        report,
        left_on="dataset2",
        right_on="Dataset",
        how="inner",
        suffixes=("_1", "_2"),
    )
    .drop(columns=["dataset1", "dataset2"])
    .sort_values("TAPE", ascending=True)
)

# %%
targets[
    ["Target label_1", "Target label_2", "Biosample name_1", "Biosample name_2"]
].head()

# %%
hepg2_targets = (
    targets[
        (targets["Biosample name_1"] == "HepG2")
        & (targets["Biosample name_2"] == "HepG2")
    ]
    .reset_index(drop=True)
    .head(200)
)

# %%
hepg2_targets[["Target label_1", "Target label_2", "blastp", "simpson"]].to_csv(
    "hepg2_targets.csv", index=False
)

# %%
data["blastp"].describe()

# %%
hepg2_targets["blastp"].describe()

# %%

data["TAPE"].describe()


# %%
hepg2_targets["TAPE"].describe()
# %%
data["simpson"].describe()

# %%
hepg2_targets["simpson"].describe()

# %%
