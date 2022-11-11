# 2.1 上位10%がなにかを出力で示す
# 上位10%をKeywordCosineとTAPEの類似度から選択する
# %%
import seaborn as sns
import matplotlib.pyplot as plt
import os
from dotenv import load_dotenv
from src.plot.util.similarity_representative import (
    get_ranked_similarity,
    get_topN_protein_pair_keyword_appended,
)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVE_DIR = os.path.join(
    PROJECT_PATH,
    "src/plot/img/select_representative_dataset_pair",
)

if not os.path.exists(SAVE_DIR):
    os.makedirs(SAVE_DIR)

# %%
# KeywordCosineとTAPEの類似度をプロット

data, _ = get_ranked_similarity()
TOPN = 100

fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(
    data=data, x="TAPE", y="Keyword", hue=data["rank"] < TOPN, alpha=0.3, ax=ax
)
fig.savefig(os.path.join(SAVE_DIR, "TAPE_Keyword_scatter.png"))

# %%
data["label"] = data["rank"] <= TOPN
splot = sns.pairplot(
    data,
    hue="label",
    vars=["TAPE", "Keyword", "MSA"],
    plot_kws={"alpha": 0.3},
)
splot.fig.savefig(os.path.join(SAVE_DIR, "TAPE_Keyword_MSA_pairplot.png"))
data.drop("label", axis=1, inplace=True)
# %%

proteins = get_topN_protein_pair_keyword_appended(100)
print(proteins.head())
proteins.to_csv(
    os.path.join(
        SAVE_DIR,
        "top100_protein_pair.csv",
    ),
    index=False,
)

# %%
proteins["Keywords1"].map(
    lambda x: [key.strip() for key in x.split(";")]
).explode().value_counts().head(20)

# %%
proteins["Keywords2"].map(
    lambda x: [key.strip() for key in x.split(";")]
).explode().value_counts().head(20)

# %%

