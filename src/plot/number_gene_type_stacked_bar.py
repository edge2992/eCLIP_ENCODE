# 実験で検出された遺伝子の種類数を遺伝子の種類で色付けして
# 積立棒グラフで描画する
# %%
import os
import pandas as pd
import sys
from dotenv import load_dotenv
import matplotlib.pyplot as plt

print(os.getcwd())


load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
sys.path.append(PROJECT_PATH)
from src.util.bedfile import load_report
from src.plot.util.process_report import count_gene

# %%


def squash_columns(df: pd.DataFrame, threshold: int) -> pd.DataFrame:
    """値が少ないカラムはotherでまとめる"""
    df.fillna(0, inplace=True)
    columns_exclude = list(gene_count_type.max().loc[lambda x: x <= threshold].index)  # type: ignore
    result = df.loc[:, ~df.columns.isin(columns_exclude)].copy()  # type: ignore
    print("squashed columns: ", columns_exclude)
    result["others"] = list(df[columns_exclude].sum(axis=1))
    return result


def split_dataframe(df: pd.DataFrame):
    """split dataframe into 2 parts"""
    n = int(df.shape[0] / 2)
    df.sort_index(inplace=True)
    return df.iloc[:n, :], df.iloc[n:, :]


# %%
report = load_report()
gene_count_type = count_gene(report[report["Biological replicates"] == "1,2"])

# %%
gene_count_type.head()
# %%
print(gene_count_type.describe())
gene_count_type.max().sort_values(ascending=False)

# %%


dfs = split_dataframe(squash_columns(gene_count_type, 20))
fig, axes = plt.subplots(2, 1, figsize=(60, 30))
fig.subplots_adjust(left=0.05, right=0.85, bottom=0.2, hspace=0.8)
for ax, data in zip(axes, dfs):
    data.plot(kind="bar", stacked=True, ax=ax, fontsize=20)
    ax.legend(fontsize=20, bbox_to_anchor=(1, 1), loc="upper left")
    ax.set_ylim(0, 6000)
fig.suptitle(
    "Number of genes per assay when converting peaks to genes using replicate1,2",
    fontsize=30,
)
fig.savefig(
    os.path.join(PROJECT_PATH, "src/plot/img", "number_gene_type_stacked_bar.png")
)

# %%
