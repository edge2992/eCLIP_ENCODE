# %%
# 遺伝子ごとにbedfileに収録されている結合部位の数を取得してプロットする
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm

PROJECT_PATH = "/mnt/H/MYWORK/eCLIP_ENCODE"
sys.path.append(PROJECT_PATH)
from src.plot.util.bedfile import get_file_path, count_file_length

# %%
def calc_all_binding_count(report: pd.DataFrame):
    """全ての遺伝子の結合部位の数を取得する"""
    """ dict gene -> -> biosample -> accession"""
    genes = report["Target label"].unique()
    d = {}
    for gene in tqdm(genes):
        report_ex = report[report["Target label"] == gene]
        biosamples = report_ex["Biosample name"].unique()
        dfs = {}
        for sample in biosamples:
            report_ex_sample = report_ex[report_ex["Biosample name"] == sample]
            dd = {}
            for _, row in report_ex_sample.iterrows():
                # label = "{} replicate_{}\n{}".format(
                #     row["Biosample name"],
                #     row["Biological replicates"],
                #     row["Accession"],
                # )
                label = row["Accession"]
                dd[label] = count_file_length(get_file_path(row))
            dfs[sample] = dd
        d[gene] = dfs
    return d


def format_N_binding_per_experiment(report: pd.DataFrame):
    """実験ごとに結合部位の数を整形する"""
    result = calc_all_binding_count(report)
    dd = {}
    for dataset in report["Dataset"].unique():
        sample_ex = report[report["Dataset"] == dataset]
        assert len(sample_ex["Target label"].unique()) == 1
        gene = sample_ex["Target label"].unique()[0]
        d = {}
        for _, row in sample_ex.iterrows():
            label = "replicates_{}".format(row["Biological replicates"])
            d[label] = result[gene][row["Biosample name"]][row["Accession"]]
        dd[dataset] = d
    return pd.DataFrame(dd).T


# %%
# report.tsvを読み取って、ファイルの対応関係を作る
report = pd.read_table(os.path.join(PROJECT_PATH, "data", "report.tsv"), skiprows=1)
report.sort_values("Biological replicates", inplace=True)


# %%
# 実験ごとに結合部位の数を数える
count_binding = format_N_binding_per_experiment(report)

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
