# 遺伝子ごとに相互作用が検出された実験の種類数を調べる
# replicate1, 2で見る
# あとから別のものが追加された時に対応できるようにコードを書く
# %%
import os
from typing import Dict, List
import pandas as pd
import sys
from tqdm import tqdm
from joblib import Parallel, delayed
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


def get_unique_gene_list(row: pd.Series) -> List[str]:
    """annotated bedを読み込んで、遺伝子のリストを作成する"""
    df = read_annotated_bed(get_formatted_file_path(row, FormatStrategy.MAX))
    return df["gene"].dropna().unique().tolist()  # type: ignore


def create_accession_gene_dict(report: pd.DataFrame):
    # Dict[accession, List[gene]]を作成する
    def do_parallel(row, how):
        df = read_annotated_bed(get_formatted_file_path(row, how))
        label = row["Dataset"]
        return label, df["gene_id"].dropna().unique().tolist()  # type: ignore

    # 並列化
    accession_gene = Parallel(n_jobs=5, verbose=3)(
        delayed(do_parallel)(row, FormatStrategy.MAX) for _, row in report.iterrows()
    )

    if accession_gene is None:
        raise ValueError("accession_gene is None")

    # format to dict
    return {key: value for key, value in accession_gene}


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


def count_experiments_by_gene(report: pd.DataFrame):
    """遺伝子ごとに検出された実験数をカウントする"""
    accession_gene = create_accession_gene_dict(report)
    gene_accession = create_gene_accession_dict(accession_gene)

    # count
    gene_count = pd.Series(
        {gene: len(accession_list) for gene, accession_list in gene_accession.items()}
    ).sort_values(ascending=False)
    return gene_count


# %%
report = load_report()

result = []
for replicate in report["Biological replicates"].unique():
    data = count_experiments_by_gene(
        report[report["Biological replicates"] == replicate]  # type: ignore
    )
    data.name = replicate
    result.append(data)

# %%
data = pd.DataFrame(result).fillna(0).T.add_prefix("replicate_")

PLOT_DIR = os.path.join(
    PROJECT_PATH,
    "src",
    "plot",
    "img",
    "number_interacting_protein_by_gene_distplot",
)
if not os.path.exists(PLOT_DIR):
    os.makedirs(PLOT_DIR)

# %%
# 横に並べるなら軸を揃えたい
replicate_num = data.shape[1]

fig, axes = plt.subplots(1, replicate_num, figsize=(7 * replicate_num, 7))
fig.subplots_adjust(bottom=0.2)

for ax, (label, values) in zip(axes, data.items()):
    sns.histplot(values[values > 0], ax=ax)
    ax.set_xlabel("# experiments in which the gene is detected", fontsize=14)
    ax.set_ylabel("Counts of gene", fontsize=14)
    ax.set_title(label, fontsize=16)  # type: ignore
    ax.set_xlim(0, 250)
    ax.set_ylim(0, 12000)

fig.suptitle(
    "Distribution of the number of eCLIP experiments in which the gene is detected",
    fontsize=18,
    y=0.1,
)

fig.savefig(os.path.join(PLOT_DIR, "compare_distplot.png"))

# %%
# replicate1,2だけを拡大して描画する
fig, ax = plt.subplots(figsize=(8, 5))
fig.subplots_adjust(bottom=0.3)

sns.histplot(
    data["replicate_1,2"][data["replicate_1,2"] > 0], ax=ax, label="replicate_1,2"
)
ax.set_xlabel("# experiments in which the gene is detected", fontsize=12)
ax.set_ylabel("Counts of gene", fontsize=12)
ax.set_title("replicate_1,2", fontsize=14)  # type: ignore

fig.suptitle(
    "Distribution of the number of eCLIP experiments\n in which the gene is detected",
    fontsize=14,
    y=0.1,
)
fig.savefig(os.path.join(PLOT_DIR, "replicate12_distplot.png"))

# %%
