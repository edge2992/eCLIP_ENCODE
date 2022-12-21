# stringdbを使って、jaccardとsimpsonの相互作用スコアが高いものを確認する
# Jaccardスコアが高い、低いの定義を作ってケースごとに可視化する

# %%
import os

import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.plot.interaction_metrics.representative import target_report
from src.util.metrics import (
    Condition,
    ConditionAnd,
    ConditionGt,
    ConditionLt,
    Metrics,
)
from src.util.similarity_strategy import DirectStringScore, Jaccard, Simpson
from src.util.metrics.keyword import dataframe_fisher_exact

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

# Jaccardスコアが高い、低いの定義を作ってケースごとに可視化する
# %%

FILEPREFIX = f"{BIOSAMPLE}_threshold_gene_num_{THRESHOLD_GENE_NUM}"


def save_fisher_exact_test(condition: Condition, data: pd.DataFrame):
    sample = dataframe_fisher_exact(condition, data)
    sample.to_csv(
        os.path.join(
            SAVEDIR,
            FILEPREFIX + str(condition).replace(" ", "_") + ".csv",
        ),
    )
    return sample


# %%
# HIGH_JACCARD_THRESHOLD = 0.2
# LOW_JACCARD_THRESHOLD = 0.1
# HIGH_STRING_DB_THRESHOLD = 0.9
# LOW_STRING_DB_THRESHOLD = 0.1
high_stringdb_score = ConditionGt("stringdb_score", 0.9)
high_jaccard = ConditionGt("Jaccard", 0.2)
low_stringdb_score = ConditionLt("stringdb_score", 0.1)
low_jaccard = ConditionLt("Jaccard", 0.1)

#
for condition in [
    ConditionAnd([high_stringdb_score, high_jaccard]),
    ConditionAnd([high_stringdb_score, low_jaccard]),
    ConditionAnd([low_stringdb_score, high_jaccard]),
    ConditionAnd([low_stringdb_score, low_jaccard]),
]:
    print(condition, f"#{data[condition(data)].shape[0]}")

# %%
# Case 1: stringdb_scoreが高くて、jaccardスコアが高いものの特徴はどうなっているか
# %%
# Case 2: stringdb_scoreが高くて、jaccardスコアが低いものの特徴はどうなっているか
save_fisher_exact_test(ConditionAnd([high_stringdb_score, low_jaccard]), data)
# %%
# Case 3: stringdb_scoreが低くてJaccardスコアが高いものの特徴はどうなっっているか
save_fisher_exact_test(ConditionAnd([low_stringdb_score, high_jaccard]), data)

# %%
# Case 4: stringdb_scoreが低くてJaccardスコアが低いものの特徴はどうなっているか
save_fisher_exact_test(ConditionAnd([low_stringdb_score, low_jaccard]), data)

# %%
