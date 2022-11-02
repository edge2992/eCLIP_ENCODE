# 共通して相互作用する遺伝子数をヒートマップで表示する
# replicate1, 2を使用
# 遺伝子数が1000件以上のアッセイ で行う

# %%
import os
import pandas as pd
from typing import Callable, List, Dict
import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv
from src.plot.util.process_report import (
    count_gene,
    get_gene_ids,
    label_protein_biosample,
)
from src.util.bedfile import load_report
from tqdm import tqdm

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

# %%

report = load_report()


def assay_with_many_gene(
    report: pd.DataFrame,
    threshold: int,
    prepare_label: Callable[[pd.DataFrame], pd.Series] = lambda x: x["Accession"],
) -> List[str]:
    """遺伝子数が多いアッセイを取得する"""
    data: pd.Series[int] = count_gene(report, prepare_label).sum(axis=1)
    intended_indexes: List[str] = sorted(list(data[data > threshold].index))
    print("{} -> {}".format(data.shape[0], len(intended_indexes)))
    return intended_indexes


def count_interection(accession_genes: Dict[str, List[str]]) -> pd.DataFrame:
    columns = sorted(list(accession_genes.keys()))
    data = pd.DataFrame(columns=columns, index=columns, dtype=int)
    for i, (k1, v1) in tqdm(enumerate(accession_genes.items())):
        for j, (k2, v2) in enumerate(accession_genes.items()):
            data.loc[k1, k2] = len(set(v1) & set(v2))
    return data


# %%
# かぶりを数える
accession_genes: Dict[str, List[str]] = get_gene_ids(
    report[report["Biological replicates"] == "1,2"], label_protein_biosample
)

count_all = count_interection(accession_genes)
# %%
# 抽出する
intended_labels = assay_with_many_gene(
    report[report["Biological replicates"] == "1,2"], 1000, label_protein_biosample
)

data = count_all.loc[intended_labels, intended_labels]

fig, ax = plt.subplots(figsize=(50, 50))
sns.heatmap(
    data.astype(int), ax=ax, cmap="Blues", annot=True, fmt="d", linewidth=0.5, vmax=1000
)
fig.savefig(
    os.path.join(
        PROJECT_PATH, "src", "plot", "img", "number_commonly_gene_heatmap.png"
    ),
)
plt.clf()
plt.close()

# %%
# fig, ax = plt.subplots(figsize=(200, 200))
# sns.heatmap(
#     count_all.astype(int),
#     ax=ax,
#     cmap="Blues",
#     annot=True,
#     fmt="d",
#     linewidth=0.5,
#     vmax=1000,
# )
# fig.savefig(
#     os.path.join(
#         PROJECT_PATH, "src", "plot", "img", "number_commonly_gene_heatmap_all.png"
#     ),
# )
# plt.clf()
# plt.close()
