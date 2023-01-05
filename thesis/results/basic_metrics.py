# %% import packages
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.eclip.dataset import Dataset
from src.util.bedfile import load_replicateIDR_report
from thesis.utils.matplotlib_format import MATPLOTLIB_CONFIG

load_dotenv()

THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
SAVEDIR = os.path.join(THESIS_FIG_PATH, "results")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

for key, value in MATPLOTLIB_CONFIG.items():
    plt.rcParams[key] = value

# %% peak_gene_scatterplot
# 遺伝子数と結合部位数の散布図を描画する
report = load_replicateIDR_report(drop_duplicates=True)
d = []
for index, row in report.iterrows():
    dataset = Dataset(row)
    d.append(
        {
            "dataset": dataset.dataset,
            "gene": len(dataset.genes),
            "peak": dataset.peaks_len,
            "biosample": dataset.biosample,
        }
    )
data = pd.DataFrame(d)

# %% タンパク質の件数
report["Target label"].nunique()
# %% plot

fig, ax = plt.subplots(1, 1)
sns.scatterplot(
    data, x="gene", y="peak", hue="biosample", ax=ax, s=80, edgecolor="none"
)
ax.set_xlabel("number of unique genes")
ax.set_ylabel("number of unique peaks")
ax.set_title(
    f"Pearson's correration coefficient: {data['gene'].corr(data['peak']):.3f}"
)
ax.grid(which="major", ls="-")
fig.savefig(os.path.join(SAVEDIR, "peak_gene_scatterplot.pdf"), dpi=300)

# %%
