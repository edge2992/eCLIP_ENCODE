# tapeをTSNEで次元圧縮してプロットする

# %%
from typing import List, Union
import numpy as np
import pandas as pd
import os
from dotenv import load_dotenv

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
EMBEDDING_PATH = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_tape.npz"

# https://builtin.com/data-science/tsne-python
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib.pyplot as plt
import seaborn as sns

# %%


def load_embedding():
    arrays = np.load(EMBEDDING_PATH, allow_pickle=True)
    return {key: value[()]["avg"] for key, value in arrays.items()}


# %%
# タンパク質ごとにラベルを付ける

KEYWORD_FILE = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/keyword_list.tsv"
KEYWORD_PROTEIN = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_keyword.tsv"


def label_protein(drop_list: list = []) -> pd.DataFrame:
    """タンパク質のラベルの変換テーブルを作成する"""
    keyword_explode: pd.Series = (
        pd.read_table(KEYWORD_PROTEIN)["Keyword ID"]
        .map(lambda x: [key.strip() for key in x.split(";")])
        .explode()
    )

    keyword_explode_annot = pd.merge(
        left=pd.concat(
            [
                keyword_explode,
                pd.read_table(KEYWORD_PROTEIN)[["From", "Entry", "Entry Name"]],
            ],
            axis=1,
        ),
        right=pd.read_table(KEYWORD_FILE)[["Keyword ID", "Name", "Category"]],
        on="Keyword ID",
        how="left",
    )
    bio_process_result = keyword_explode_annot[
        keyword_explode_annot["Category"] == "Biological process"
    ]
    left = bio_process_result[~bio_process_result["Name"].isin(drop_list)]
    right = (
        left["Name"]
        .value_counts()
        .reset_index()
        .rename(columns={"index": "Name", "Name": "Count"})
    )
    result = (
        left.join(right.set_index("Name"), on="Name")
        .sort_values("Count", ascending=False)
        .groupby("From")
        .first()
        .reset_index()
    )
    print(
        "Protein count: {} -> {}".format(
            pd.read_table(KEYWORD_PROTEIN)["From"].nunique(), result["From"].nunique()
        )
    )
    return result


def convert_label(
    proteins: List[str], drop_list: List[str], threshold: Union[None, int] = None
) -> List[str]:
    """proteinのリストをラベルに変換する"""
    label_table = label_protein(drop_list)
    if threshold is not None:
        vv = label_table["Name"].value_counts()
        label_table = label_table[label_table["Name"].isin(vv[vv >= threshold].index)]
    print("label count: {}".format(label_table["Name"].nunique()))
    print(label_table["Name"].value_counts())
    name_dict = {
        "sp|" + row["Entry"] + "|" + row["Entry Name"]: row["Name"]
        for _, row in label_table.iterrows()
    }
    return [name_dict[key] if key in name_dict.keys() else "other" for key in proteins]


# %%

DROP_LIST = [
    "mRNA processing",
    "Transcription",
    "Transport",
    "mRNA Transport" "Innate immunity",
]

# %%
arrays = load_embedding()

# %%

pca = PCA(n_components=3)
pca_result = pca.fit_transform(list(arrays.values()))
labels = convert_label(list(arrays.keys()), DROP_LIST, threshold=5)

fig, ax = plt.subplots(figsize=(16, 10))
sns.scatterplot(x=pca_result[:, 0], y=pca_result[:, 1], ax=ax, hue=labels)
ax.set_xlabel("PC1", fontsize=12)
ax.set_ylabel("PC2", fontsize=12)
fig.savefig(os.path.join(PROJECT_PATH, "src/plot/img", "tape_tsne", "pca_scatter.png"))

# %%
# plot tsne

tsne = TSNE(n_components=2, verbose=1, perplexity=40, n_iter=300)
tsne_results = tsne.fit_transform(np.array(list(arrays.values())))
labels = convert_label(list(arrays.keys()), DROP_LIST, threshold=5)

# %%
fig, ax = plt.subplots(figsize=(16, 10))
sns.scatterplot(x=tsne_results[:, 0], y=tsne_results[:, 1], ax=ax, hue=labels)
ax.set_xlabel("tsne-1", fontsize=12)
ax.set_ylabel("tsne-2", fontsize=12)
fig.savefig(os.path.join(PROJECT_PATH, "src/plot/img", "tape_tsne", "tsne_scatter.png"))

# %%
tsne_df = pd.DataFrame(
    {"tsne-1": tsne_results[:, 0], "tsne-2": tsne_results[:, 1], "label": labels}
)
fig, ax = plt.subplots(figsize=(16, 10))
sns.scatterplot(
    tsne_df[tsne_df["label"] != "other"],
    x="tsne-1",
    y="tsne-2",
    hue="label",
    ax=ax,
    palette=sns.color_palette("hls", len(tsne_df["label"].unique()) - 1),
)
ax.set_xlabel("tsne-1", fontsize=12)
ax.set_ylabel("tsne-2", fontsize=12)
fig.savefig(
    os.path.join(
        PROJECT_PATH, "src/plot/img", "tape_tsne", "tsne_scatter_exclude_other.png"
    )
)

# %%
