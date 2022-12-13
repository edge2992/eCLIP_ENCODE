# keywordごとにstringdbのスコアの分布を見てみる

# %%
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.plot.interaction_metrics.representative import (
    convert_to_dict_exp_pair_by_keyword,
    target_report,
    metrics,
)
from src.util.similarity_strategy import Jaccard, Simpson, TAPE, DirectStringScore

sns.set(font_scale=1.4)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(
    PROJECT_PATH, "src/plot/img", "foldseek_stringdb_investigation", "keyword_stringdb"
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

    # %%

THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
interaction_strategies = {
    "simpson": lambda report: Simpson(report),
    "jaccard": lambda report: Jaccard(report),
}
protein_strategies = {
    "TAPE": lambda report: TAPE(report),
    "stringdb_score": lambda report: DirectStringScore(report, metrics="score"),
    "stringdb_ascore": lambda report: DirectStringScore(report, metrics="ascore"),
    "stringdb_escore": lambda report: DirectStringScore(report, metrics="escore"),
    "stringdb_tscore": lambda report: DirectStringScore(report, metrics="tscore"),
}

data = metrics(report, protein_strategies, interaction_strategies)

# %%
keyword_experiment_pair = convert_to_dict_exp_pair_by_keyword(data)

# %%
data.head()
plot_data = pd.DataFrame()
for keyword, value in keyword_experiment_pair.items():
    print(keyword, len(value))
    print(data.iloc[value, :].head())
    sample = (
        data.iloc[value, :][
            ["stringdb_score", "stringdb_ascore", "stringdb_escore", "stringdb_tscore"]
        ]
        .copy()
        .reset_index()
    )
    sample["label"] = keyword
    plot_data = pd.concat([plot_data, sample], axis=0)

# %%
print(plot_data.shape)
print(plot_data.head())
fig, ax = plt.subplots(figsize=(100, 10))
sns.boxplot(plot_data, x="label", y="stringdb_score", ax=ax)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
fig.tight_layout()
fig.savefig(os.path.join(SAVEDIR, "boxplot_stringdb_score.png"))

# TODO: 全てと比較する
# labelの文字を大きくする 並び替えをする
# 色を付けない
# ラベルにデータ数を含める


# %%
