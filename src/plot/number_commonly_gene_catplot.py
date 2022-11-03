# タンパク質の種類ごとに
# 他のアッセイとの遺伝子の被りの個数を箱ひげ図で表示する

# %%
import os
import seaborn as sns
import matplotlib.pyplot as plt
from dotenv import load_dotenv
from src.plot.util.split import split_dataframe
from src.plot.util.process_intersect_gene import count_interection
from src.plot.util.process_report import label_protein_biosample
from src.util.bedfile import load_report

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

# %%


# %%
# 被りを数える

report = load_report()

count_all = count_interection(
    report[report["Biological replicates"] == "1,2"], label_protein_biosample
)


# %%
# 全てのタンパク質でそれぞれ箱ひげ図を書く
SPLIT_NUM = 10

fig, axes = plt.subplots(SPLIT_NUM, 1, figsize=((400 / SPLIT_NUM), 20 * SPLIT_NUM))
dfs = split_dataframe(count_all, SPLIT_NUM)

color_dict = {
    target: sns.color_palette()[i]  # type: ignore
    for i, target in enumerate(count_all.columns.map(lambda x: x.split()[1]).unique())
}

for ax, data in zip(axes, dfs):
    palette = [color_dict[label.split()[1]] for label in data.index]
    sns.boxplot(data=data.T, ax=ax, palette=palette)

fig.savefig(
    os.path.join(
        PROJECT_PATH, "src", "plot", "img", "number_commonly_gene_violinplot.png"
    )
)

# %%
# 箱ひげ図を一つ書く
target_protein = "BUD13"

report_ex = report[
    (report["Target label"] == target_protein)
    & (report["Biological replicates"] == "1,2")
]
n_catplot = report_ex.shape[0]  # seriesでも大丈夫?
data = count_all.loc[label_protein_biosample(report_ex)]
sns.catplot(data=data.T, kind="violin")

# %%
