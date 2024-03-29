# %%
# タンパク質 (RBP) ごとにgene annotated bedfileに収録されている結合部位の数を取得してプロットする
import pandas as pd
import os
import sys
import seaborn as sns
from dotenv import load_dotenv


load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
sys.path.append(PROJECT_PATH)
from src.util.bed_format_strategy import FormatStrategy
from src.util.bedfile import count_gene_nunique
from src.util.get_bed_path import get_formatted_file_path
from src.plot.util.process_by_accession import (
    convert_gene_biosample_accession_count_dict_to_df,
    create_gene_biosample_accession_count_dict,
)

# %%
# report.tsvを読み取って、ファイルの対応関係を作る
report = pd.read_table(os.path.join(PROJECT_PATH, "data", "report.tsv"), skiprows=1)
report.sort_values("Biological replicates", inplace=True)


# %%
# 実験ごとに結合部位の数を数える
count_gene = convert_gene_biosample_accession_count_dict_to_df(
    create_gene_biosample_accession_count_dict(
        report,
        lambda row: count_gene_nunique(
            get_formatted_file_path(row, FormatStrategy.MAX)  # type: ignore
        ),
    ),
    report,
)
# %%
count_gene.head()

# %%
# 描画
# 細胞株（実験前提）ごとに色分けしたい
data = count_gene.join(
    report[["Dataset", "Biosample name"]].drop_duplicates().set_index("Dataset")
)
splot = sns.pairplot(data, hue="Biosample name")
splot.fig.suptitle(
    "Number of gene per experimental condition\n (generated by bedtools intersects)",
)
splot.fig.subplots_adjust(top=0.9)
splot.fig.savefig(
    os.path.join(PROJECT_PATH, "src", "plot", "img", "number_gene_pairplot.png"),
)
