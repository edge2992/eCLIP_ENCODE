# uniprotのキーワードを見て特徴を掴みたい
# Dict[keyword, List[Protein]]を作成する
# 沢山あるキーワードから、相互作用の数との関連を調べる

# %%
from typing import Union
import pandas as pd
import os
from dotenv import load_dotenv
from itertools import combinations
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns

from src.util.uniprot import keyword_count

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
KEYWORD_FILE = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/keyword_list.tsv"
KEYWORD_PROTEIN = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_keyword.tsv"

# %%
df = keyword_count()
print(df.shape)
print(df.head())
df.to_csv(os.path.join(PROJECT_PATH, "data", "keyword_count.tsv"), sep="\t")
# %%
print("---", "Domain", "---")
print(df[df["Category"] == "Domain"])
# %%
print("---", "Biological process", "---")
print(df[df["Category"] == "Biological process"])

# %%


def count_matrix(threshold: int = 0, category: Union[str, None] = None) -> pd.DataFrame:
    """keywordの組み合わせでheatmap用のデータを作成する"""
    tmp: pd.DataFrame = (
        pd.read_table(KEYWORD_PROTEIN)["Keyword ID"]
        .map(lambda x: [key.strip() for key in x.split(";")])
        .explode()
        .reset_index()
    )
    value_counts = tmp["Keyword ID"].value_counts()
    keys = value_counts[value_counts > threshold].index
    if category is not None:
        keytable = pd.read_table(KEYWORD_FILE)
        keytable = keytable[keytable["Category"] == category]
        keys = [key for key in keys if key in keytable["Keyword ID"].values]
    print("keys: ", keys[:5])

    # ２次元配列にする
    data = pd.DataFrame(columns=keys, index=keys, dtype=int, data=0)
    for _, df in tqdm(tmp.groupby("index")):
        keywords = [keyword for keyword in df["Keyword ID"].values if keyword in keys]
        for key in keywords:
            data.loc[key, key] += 1
        for a, b in combinations(keywords, 2):
            data.loc[a, b] += 1
            data.loc[b, a] += 1

    # convert columns index
    convert_dict = (
        pd.read_table(KEYWORD_FILE)[["Keyword ID", "Name"]]
        .set_index("Keyword ID", drop=True)["Name"]
        .to_dict()
    )
    data.rename(convert_dict, axis=0, inplace=True)
    data.rename(convert_dict, axis=1, inplace=True)
    return data


# %% 作図

fig, ax = plt.subplots(figsize=(20, 20))
sns.heatmap(
    count_matrix(10), ax=ax, cmap="Blues", annot=True, fmt="d", linewidth=0.5, vmax=50
)
fig.suptitle("Keyword count matrix", fontsize=20)
fig.savefig(
    os.path.join(PROJECT_PATH, "src/plot/img", "mostcommon_keywords", "heatmap.png")
)

# %%

fig, ax = plt.subplots(figsize=(20, 20))
sns.heatmap(
    count_matrix(0, "Domain"),
    ax=ax,
    cmap="Blues",
    annot=True,
    fmt="d",
    linewidth=0.5,
    annot_kws={"size": 16},
)
ax.tick_params(axis="x", labelsize=16)
ax.tick_params(axis="y", labelsize=16)
fig.suptitle("Domain", fontsize=20)
fig.tight_layout()
fig.savefig(
    os.path.join(
        PROJECT_PATH, "src/plot/img", "mostcommon_keywords", "heatmap_domain.png"
    )
)

# %%

fig, ax = plt.subplots(figsize=(20, 20))
sns.heatmap(
    count_matrix(5, "Biological process"),
    ax=ax,
    cmap="Blues",
    annot=True,
    fmt="d",
    linewidth=0.5,
    annot_kws={"size": 16},
)
ax.tick_params(axis="x", labelsize=16)
ax.tick_params(axis="y", labelsize=16)
fig.suptitle("Biological process", fontsize=20)
fig.tight_layout()
fig.savefig(
    os.path.join(
        PROJECT_PATH,
        "src/plot/img",
        "mostcommon_keywords",
        "heatmap_biological_process.png",
    )
)

# %%
