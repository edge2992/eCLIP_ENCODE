# %%
# eCLIPの実験ごとの相互作用遺伝子と結合部位の種類数の散布図。
# bed+narrowPeakファイルをENDSEMBLのGFFファイルとインターセクトさせて
# それぞれの結合部位を遺伝子に変換した。
import pandas as pd
import os
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from dotenv import load_dotenv


load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
sys.path.append(PROJECT_PATH)
from src.util.bed_format_strategy import FormatStrategy

from src.util.bedfile import count_file_length, count_gene_nunique
from src.util.get_bed_path import get_file_path, get_formatted_file_path
from src.plot.util.process_by_accession import (
    convert_gene_biosample_accession_count_dict_to_df,
    create_gene_biosample_accession_count_dict,
)

# %%
# report.tsvを読み取って、ファイルの対応関係を作る
report = pd.read_table(os.path.join(PROJECT_PATH, "data", "report.tsv"), skiprows=1)
report.sort_values("Biological replicates", inplace=True)

# %%
# 実験ごとに遺伝子の数を数える
count_gene = convert_gene_biosample_accession_count_dict_to_df(
    create_gene_biosample_accession_count_dict(
        report,
        lambda row: count_gene_nunique(
            get_formatted_file_path(row, FormatStrategy.MAX)
        ),
    ),
    report,
)
# %%
# 実験ごとに結合部位の数を数える
count_binding = convert_gene_biosample_accession_count_dict_to_df(
    create_gene_biosample_accession_count_dict(
        report, lambda row: count_file_length(get_file_path(row))
    ),
    report,
)
# %%
cg = count_gene.add_prefix("gene_")
cb = count_binding.add_prefix("peak_")
count_gene_binding = pd.concat([cg, cb], axis=1)
print(count_gene_binding.head())

# %%
# 描画
# 細胞株（実験前提）ごとに色分けしたい

data = count_gene_binding.join(
    report[["Dataset", "Biosample name"]].drop_duplicates().set_index("Dataset")
)
# %%

fig, axes = plt.subplots(1, 3, figsize=(20, 7))

for i, (ax, gene_col, bindingsite_col) in enumerate(zip(axes, cg.columns, cb.columns)):
    sns.scatterplot(data, x=bindingsite_col, y=gene_col, hue="Biosample name", ax=ax)
    ax.legend(loc="lower right")
    ax.set_xlabel(bindingsite_col, fontsize=16)
    ax.set_ylabel(gene_col, fontsize=16)
    ax.set_title(f"correlation: {data[gene_col].corr(data[bindingsite_col]):.3f}")

fig.subplots_adjust(bottom=0.2)
fig.suptitle(
    "Scatter plot of the number of types of interacting genes and binding sites for each eCLIP experiment.\n"
    + "The bed+narrowPeak file was intersected with the ENSEMBL GFF file to convert each binding site into a gene.",
    fontsize=16,
    y=0.1,
)
fig.savefig(
    os.path.join(
        PROJECT_PATH, "src", "plot", "img", "number_bindingsite_gene_pairplot.png"
    ),
)

# %%
