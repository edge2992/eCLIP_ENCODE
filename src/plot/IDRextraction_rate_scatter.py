# IDR によるピークの抽出率とIDR抽出後のピークの件数を散布図で可視化する
# IDRの抽出率が高いもので、タンパク質の種類のラベルを確認する (散布図で特定のラベルを色分け表示)

# 減少率の計算式
# Datasetごとに計算する
# Simpson係数に類似
# | # peak in replicate1, 2 | / min(| # peak in replicate1 |, | # peak in replicate2 |)


# %%

from src.plot.util.process_report import count_binding, count_gene
from src.util.bedfile import load_report, load_replicateIDR_report
from src.util.uniprot import keyword_count
from src.util.uniprot import load_splited_keyword_report

import pandas as pd
import os
from typing import List
from dotenv import load_dotenv
import seaborn as sns
import matplotlib.pyplot as plt


load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]

# %%


def IDR_extraction_rate(group: pd.DataFrame):
    return {
        "IDR_extraction_rate": group.min() / group.max(),
        "#peak IDR extracted": group.min(),
    }


number_peak_all = count_binding(load_report())
number_peak_all.name = "#peak"
# number_peak = count_binding(load_replicateIDR_report(), lambda df: df["Dataset"])
number_gene = count_gene(load_replicateIDR_report(), lambda df: df["Dataset"])

# %%
extraction_rate = (
    load_report()
    .join(number_peak_all, on="Accession")
    .groupby("Dataset")["#peak"]
    .apply(IDR_extraction_rate)
    .unstack()
)

data = pd.concat([extraction_rate, number_gene.rename("#gene IDR extracted")], axis=1)

data.head()

# %%

sns.set(font_scale=1.4)

fig, axes = plt.subplots(1, 3, figsize=(20, 7))

sns.scatterplot(
    data,
    x="#peak IDR extracted",
    y="#gene IDR extracted",
    ax=axes[0],
    hue="IDR_extraction_rate",
)
sns.scatterplot(
    data,
    x="#peak IDR extracted",
    y="IDR_extraction_rate",
    ax=axes[1],
    hue="#gene IDR extracted",
)
sns.scatterplot(
    data,
    x="#gene IDR extracted",
    y="IDR_extraction_rate",
    ax=axes[2],
    hue="#peak IDR extracted",
)
fig.subplots_adjust(bottom=0.2)
fig.suptitle(
    "IDR extraction rate and number of peaks and genes extracted by IDR",
)
fig.tight_layout()
fig.savefig(
    os.path.join(
        PROJECT_PATH, "src/plot/img", "IDRextraction_rate_scatter", "raw_scatter.png"
    )
)


# %%
load_replicateIDR_report()[["Dataset", "Accession", "Biosample name", "Target label"]]

# %%


def make_hue_list(datasets: List[str], keyword_id: str) -> List[bool]:
    """keyword_idを持つタンパク質を実験しているデータセットをTrueにするリストを作成する"""
    keywords = load_splited_keyword_report()

    protein_list = (
        keywords[keywords["Keyword ID"] == keyword_id]["From"].unique().tolist()
    )
    report = load_report()
    intended_datasets = report[report["Target label"].isin(protein_list)][
        "Dataset"
    ].to_list()
    return [dataset in intended_datasets for dataset in datasets]


def plot_scatter(data: pd.DataFrame, keyword_id: str, keyword_name: str):
    hue_list = [
        keyword_name if x else "other"
        for x in make_hue_list(list(data.index), keyword_id)
    ]
    palette = {keyword_name: sns.color_palette()[3], "other": sns.color_palette()[0]}  # type: ignore

    fig, axes = plt.subplots(1, 3, figsize=(20, 7))

    sns.scatterplot(
        data,
        x="#peak IDR extracted",
        y="#gene IDR extracted",
        ax=axes[0],
        palette=palette,
        hue=hue_list,
    )
    axes[0].legend(loc="lower right")
    sns.scatterplot(
        data,
        x="#peak IDR extracted",
        y="IDR_extraction_rate",
        ax=axes[1],
        palette=palette,
        hue=hue_list,
    )
    axes[1].legend(loc="lower right")
    sns.scatterplot(
        data,
        x="#gene IDR extracted",
        y="IDR_extraction_rate",
        ax=axes[2],
        palette=palette,
        hue=hue_list,
    )
    axes[2].legend(loc="lower right")
    fig.subplots_adjust(bottom=0.2)
    fig.suptitle(
        "IDR extraction rate and number of peaks and genes "
        + "extracted by IDR [{},  #{}]".format(
            keyword_name, sum([x == keyword_name for x in hue_list])
        )
    )
    fig.tight_layout()
    return fig


# %%


for i, (index, row) in enumerate(keyword_count().head(38).iterrows()):
    fig = plot_scatter(data, str(index), row["Name"])
    fig.savefig(
        os.path.join(
            PROJECT_PATH,
            "src/plot/img",
            "IDRextraction_rate_scatter",
            "{}_{}_scatter.png".format(i, "_".join(str(row["Name"]).split())),
        )
    )
    plt.clf()
    plt.close()

# %%
