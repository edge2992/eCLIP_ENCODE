# %%
# タンパク質 (RBP) ごとにbedfileに収録されている結合部位の数を取得してプロットする
import pandas as pd
import os
import sys
import seaborn as sns
from dotenv import load_dotenv

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]
sys.path.append(PROJECT_PATH)
from src.util.bedfile import count_file_length
from src.util.get_bed_path import get_file_path
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

count_binding = convert_gene_biosample_accession_count_dict_to_df(
    create_gene_biosample_accession_count_dict(
        report, lambda row: count_file_length(get_file_path(row))
    ),
    report,
)

# %%
# 描画
# 細胞株（実験前提）ごとに色分けしたい
data = count_binding.join(
    report[["Dataset", "Biosample name"]].drop_duplicates().set_index("Dataset")
)
splot = sns.pairplot(data, hue="Biosample name")
splot.fig.suptitle(
    "Number of binding sites per experimental condition\n (number of lines in bedfile)",
)
splot.fig.subplots_adjust(top=0.9)
splot.fig.savefig(
    os.path.join(PROJECT_PATH, "src", "plot", "img", "number_bindingsite_pairplot.png"),
)

# %%
