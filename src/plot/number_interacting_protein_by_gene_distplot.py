# 遺伝子ごとに相互作用が検出された実験の種類数を調べる
# replicate1, 2で見る
# あとから別のものが追加された時に対応できるようにコードを書く
# %%
import os
from typing import Dict, List
import pandas as pd
import sys
from tqdm import tqdm
from dotenv import load_dotenv
import seaborn as sns
import matplotlib.pyplot as plt

print(os.getcwd())


load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
sys.path.append(PROJECT_PATH)
from src.util.bedfile import load_report, read_annotated_bed
from src.util.bed_format_strategy import FormatStrategy
from src.util.get_bed_path import get_formatted_file_path

# %%

# %%


def get_unique_gene_list(row: pd.Series) -> List[str]:
    """annotated bedを読み込んで、遺伝子のリストを作成する"""
    df = read_annotated_bed(get_formatted_file_path(row, FormatStrategy.MAX))
    return df["gene"].dropna().unique().tolist()  # type: ignore


def create_accession_gene_dict(report: pd.DataFrame):
    # Dict[accession, List[gene]]を作成する
    result: Dict[str, List[str]] = {}
    for _, row in tqdm(report.iterrows()):
        df = read_annotated_bed(get_formatted_file_path(row, FormatStrategy.MAX))
        result[row["Dataset"]] = df["gene_id"].unique().tolist()
    return result


def create_gene_accession_dict(accession_gene: Dict[str, List[str]]):
    # Dict[gene, List[accession]]に変換する
    result: Dict[str, List[str]] = {}
    for key, gene_list in tqdm(accession_gene.items()):
        for gene in gene_list:
            if gene in result:
                result[gene].append(key)
            else:
                result[gene] = [key]
    return result


# %%
report = load_report()

gene_accession = create_gene_accession_dict(
    create_accession_gene_dict(report[report["Biological replicates"] == "1,2"])
)

# %%
data = pd.Series({k: len(v) for k, v in gene_accession.items()}).sort_values(
    ascending=False
)
print(data.head())

# %%

fig, ax = plt.subplots(figsize=(5, 5))
fig.subplots_adjust(top=0.8, left=0.2)

sns.histplot(data, ax=ax)

ax.set_xlabel("# experiments in which the gene is detected", fontsize=12)
ax.set_ylabel("Counts of gene", fontsize=12)
ax.set_title(
    "The number of eCLIP experiments\n"
    + "in which the gene was detected per gene.\n"
    + "The total number of gene types detected was {},\n".format(len(data))
    + "and plotted using replicate1,2.",
    fontsize=10,
    y=1.03,
)

fig.savefig(
    os.path.join(
        PROJECT_PATH,
        "src",
        "plot",
        "img",
        "number_interacting_protein_by_gene_distplot.png",
    ),
)

# %%
