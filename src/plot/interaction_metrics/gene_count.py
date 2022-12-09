# %%

from src.util.bedfile import load_replicateIDR_report
from src.plot.interaction_metrics.representative import get_geneset
from tqdm import tqdm
import pandas as pd
import seaborn as sns
import os
import matplotlib.pyplot as plt
from dotenv import load_dotenv

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "interaction_metrics", "gene_count"
)


def get_gene_count(biosample: str):
    """biosample (HepG2 or K562) ごとに相互作用する遺伝子の個数をまとめる"""
    if biosample not in ["HepG2", "K562"]:
        raise ValueError("biosample must be HepG2 or K562")
    report = load_replicateIDR_report()
    report = report[report["Biosample name"] == biosample].reset_index(drop=True)

    dataset_gene_dict = {}
    for index, row in tqdm(report.iterrows()):
        dataset = row["Dataset"]
        dataset_gene_dict[dataset] = get_geneset(dataset)

    return pd.merge(
        pd.DataFrame(
            {k: len(v) for k, v in dataset_gene_dict.items()}, index=["Gene count"]
        ).T,
        report[["Dataset", "Biosample name", "Target label"]],
        right_on="Dataset",
        left_index=True,
    ).sort_values("Gene count", ascending=False)


# %%

if os.path.exists(SAVEDIR) is False:
    os.makedirs(SAVEDIR)

sns.set(style="whitegrid", font_scale=1.5)

for biosample in ["HepG2", "K562"]:
    sample = get_gene_count(biosample)
    fig, ax = plt.subplots(figsize=(30, 5))
    sns.barplot(
        data=sample, y="Gene count", x="Target label", hue="Biosample name", ax=ax
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_title(f"{biosample} gene count")
    fig.savefig(os.path.join(SAVEDIR, f"{biosample}.png"), bbox_inches="tight")
# %%
