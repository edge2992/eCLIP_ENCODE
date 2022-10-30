# %%
import pandas as pd
import os
import matplotlib.pyplot as plt
from tqdm import tqdm

PROJECT_PATH = "/mnt/H/MYWORK/eCLIP_ENCODE"
# %%
def read_eCLIP_bed(filename):
    """eCLIPのbedfileにcolumn名を付けて読み込む"""
    columns = [
        "chrom",
        "chromStart",
        "chromEnd",
        "name",
        "score",
        "strand",
        "singleValue",
        "pValue",
        "qValue",
        "peak",
    ]
    df = pd.read_table(filename, names=columns)
    return df


def count_file_length(filename):
    """ファイルの行数を数える"""
    print(filename)
    with open(filename, "r") as f:
        return len(f.readlines())


def get_file_path(row: pd.Series):
    """report.tsvの行からファイルのパスを取得する"""
    return os.path.join(
        PROJECT_PATH,
        "data",
        row["Assay term name"],
        row["Target label"],
        row["Biosample name"].split()[0],  # adrenal gland, K562, HepG2
        row["Accession"] + ".bed",
    )


# report.tsvを読み取って、ファイルの対応関係を作る
report = pd.read_table(os.path.join(PROJECT_PATH, "data", "report.tsv"), skiprows=1)
report.sort_values("Biological replicates", inplace=True)

# %%
# singleValueとpValueの散布図を作成する
# plot
def plot_singleValue_pValue_scatter(gene: str, report: pd.DataFrame):
    """遺伝子ごとにsingleValueとpValueの散布図を作成する"""
    # assay_title = "eCLIP"
    report_ex = report[report["Target label"] == gene]
    biosamples = report_ex["Biosample name"].unique()
    dfs = []
    for sample in biosamples:
        report_ex_sample = report_ex[report_ex["Biosample name"] == sample]
        dd = {}
        for _, row in report_ex_sample.iterrows():
            label = "{} replicate_{}\n{}".format(
                row["Biosample name"], row["Biological replicates"], row["Accession"]
            )
            dd[label] = read_eCLIP_bed(get_file_path(row))
        dfs.append(dd)

    # plot as a scatter plot
    n_plot = sum([len(dd) for dd in dfs])
    n_row = len(biosamples)
    n_col = int(n_plot / n_row)
    assert n_plot == n_row * n_col, "{} != {}* {}".format(n_plot, n_row, n_col)
    fig, axes = plt.subplots(
        n_row,
        n_col,
        figsize=(15, 7.5 * len(biosamples)),
    )
    axes = axes.flatten()
    for i, dd in enumerate(dfs):
        for j, (label, df) in enumerate(dd.items()):
            ax = axes[i * n_col + j]
            ax.scatter(df["singleValue"], df["pValue"])
            ax.set_title(label)
            ax.set_xlabel("singleValue")
            ax.set_ylabel("pvalue (-log10)")
            ax.set_xlim(-1, 10)
            ax.set_ylim(1e-4, 1e3)
            ax.set_yscale("log")
    return fig


# %%
# 一つ描画してみる
gene = "UPF1"
fig = plot_singleValue_pValue_scatter(gene, report)

# %%
# 描画
# FXR2, HNRNPKで失敗している (全体をざっくりみたいだけなので無視)

SAVE_IMG_DIR = os.path.join(
    PROJECT_PATH, "src", "plot", "img", "singleValue_pValue_scatter"
)
if not os.path.exists(SAVE_IMG_DIR):
    os.makedirs(SAVE_IMG_DIR)

for gene in tqdm(report["Target label"].unique()):
    try:
        fig = plot_singleValue_pValue_scatter(gene, report)
        fig.savefig(os.path.join(SAVE_IMG_DIR, "{}.png".format(gene)))
        plt.clf()
        plt.close()
    except Exception as e:
        print(gene + " is failed")
        print(e)
        continue
