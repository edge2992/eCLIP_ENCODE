# stringdbを使って、jaccardとsimpsonの相互作用スコアが高いものを確認する

# %%
import os

import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.util.metrics import Metrics
from src.plot.interaction_metrics.representative import target_report
from src.util.similarity_strategy import (
    DirectStringScore,
    Jaccard,
    Simpson,
)
from src.plot.interaction_metrics.representative import describe_dataset_pair


sns.set(font_scale=1.4)

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

SAVEDIR = os.path.join(
    PROJECT_PATH,
    "src/plot/img",
    "foldseek_stringdb_investigation",
    "see_example_stringdb",
)

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


# %% prepare_plotting data

THRESHOLD_GENE_NUM = 1000
BIOSAMPLE = "HepG2"


report = target_report(THRESHOLD_GENE_NUM, BIOSAMPLE)
data: pd.DataFrame = Metrics(report)(
    [
        DirectStringScore(report, metrics="score"),
        DirectStringScore(report, metrics="ascore"),
        DirectStringScore(report, metrics="escore"),
        DirectStringScore(report, metrics="tscore"),
        Simpson(),
        Jaccard(),
    ]
)  # type: ignore

# %%
nonzero_data = data[data["stringdb_score"] != 0].copy()
nonzero_data.describe()

# %%
data[["simpson", "jaccard"]].describe()

HIGH_JACCARD_THRESHOLD = 0.2
LOW_JACCARD_THRESHOLD = 0.1
HIGH_STRING_DB_THRESHOLD = 0.9
LOW_STRING_DB_THRESHOLD = 0.1
LABELS = ["Target label_1", "Target label_2", "stringdb_score", "jaccard", "simpson"]


# Jaccardスコアが高い、低いの定義を作ってケースごとに可視化する

# Case 1: stringdb_scoreが高くて、jaccardスコアが高いものの特徴はどうなっているか
# %%
high_stringdb_high_jaccard = (
    nonzero_data[nonzero_data["stringdb_score"] > HIGH_STRING_DB_THRESHOLD]
    .sort_values("jaccard", ascending=False)
    .head(10)
)
print(high_stringdb_high_jaccard[LABELS])
for index, row in high_stringdb_high_jaccard.iterrows():
    describe_dataset_pair(row)
# Spliceosome関連っぽい

# %%
# Case 2: stringdb_scoreが高くて、jaccardスコアが低いものの特徴はどうなっているか

high_stringdb_low_jaccard = (
    nonzero_data[nonzero_data["stringdb_score"] > HIGH_STRING_DB_THRESHOLD]
    .sort_values("jaccard", ascending=True)
    .head(10)
)
print(high_stringdb_low_jaccard[LABELS])
for index, row in high_stringdb_low_jaccard.iterrows():
    describe_dataset_pair(row)

# %%

# Case 3: stringdb_scoreが低くてJaccardスコアが高いものの特徴はどうなっっているか

low_stringdb_high_jaccard = (
    nonzero_data[nonzero_data["stringdb_score"] < 0.2]
    .sort_values("jaccard", ascending=False)
    .head(10)
)
print(low_stringdb_high_jaccard.shape)
print(low_stringdb_high_jaccard[LABELS])
for index, row in low_stringdb_high_jaccard.iterrows():
    describe_dataset_pair(row)

# %%

# Case 4: stringdb_scoreが低くてJaccardスコアが低いものの特徴はどうなっているか

# Case 5: stringdb_scoreがないものはどうなっているか
nonzero_data[["simpson", "jaccard"]].describe()


# %%
