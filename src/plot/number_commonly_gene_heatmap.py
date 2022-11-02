# 共通して相互作用する遺伝子数をヒートマップで表示する
# replicate1, 2を使用
# 遺伝子数が1000件以上のアッセイ で行う

# %%
import os
import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv
from src.plot.util.process_intersect_gene import assay_with_many_gene, count_interection
from src.plot.util.process_report import (
    label_protein_biosample,
)
from src.util.bedfile import load_report

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

# %%

report = load_report()

count_all = count_interection(
    report[report["Biological replicates"] == "1,2"], label_protein_biosample
)
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

# %%
# 全体を表示 20MB
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

# %%
