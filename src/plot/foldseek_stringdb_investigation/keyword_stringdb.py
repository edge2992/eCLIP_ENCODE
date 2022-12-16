# keywordごとにstringdbのスコアの分布を見てみる

# %%
import os

import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.plot.foldseek_stringdb_investigation.plot_utils import (
    plot_boxplot_by_keyword,
)
from src.plot.interaction_metrics.representative import (
    target_report,
)
from src.util.metrics import Metrics, KeywordConfidence, ConditionGt
from src.util.similarity_strategy import TAPE, DirectStringScore, Jaccard, Simpson
from typing import List

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

data: pd.DataFrame = Metrics(report)(
    [
        TAPE(),
        DirectStringScore(report, metrics="score"),
        DirectStringScore(report, metrics="ascore"),
        DirectStringScore(report, metrics="escore"),
        DirectStringScore(report, metrics="tscore"),
        Simpson(),
        Jaccard(),
    ]
)  # type: ignore

# %%
# keyword_experiment_pair = convert_to_dict_exp_pair_by_keyword(data)
keyword_confidence = KeywordConfidence(data)
keyword_experiment_pair = keyword_confidence.keywords_dict
# %%


def plotting_data(data: pd.DataFrame, keyword_experiment_pair, metrics: List):
    """plot_boxplot_by_keywordで作図するデータを作成する

    Args:
        data (pd.DataFrame): from src.plot.interaction_metrics.representativeのmetricsで作成したデータ
        keyword_experiment_pair (_type_): convert_to_dict_exp_pair_by_keywordで作成したデータ
        metrics (List): 対象とするmetrics

    Returns:
        _type_: _description_
    """
    from src.eclip import Dataset, Compare

    plot_data = pd.DataFrame()
    tmp = data.copy()
    tmp["label"] = data.apply(
        lambda row: Compare(
            [Dataset(row["Dataset_1"]), Dataset(row["Dataset_2"])]
        ).label,
        axis=1,
    )
    tmp.set_index("label", inplace=True, drop=True)
    for keyword, value in keyword_experiment_pair.items():
        sample = tmp.loc[value, :][metrics].copy().reset_index()
        sample["label"] = f"{keyword} (#{len(value)})"
        plot_data = pd.concat([plot_data, sample], axis=0)
    return plot_data


investigation_metrics = [
    "stringdb_score",
    "stringdb_ascore",
    "stringdb_escore",
    "stringdb_tscore",
]
plot_data = plotting_data(data, keyword_experiment_pair, investigation_metrics)

for met in investigation_metrics:
    fig = plot_boxplot_by_keyword(plot_data, met)
    fig.savefig(os.path.join(SAVEDIR, f"boxplot_{met}.png"), bbox_inches="tight")

# %%
from tqdm import tqdm

condition = ConditionGt("stringdb_score", 0.5)
dd = []
for key in tqdm(keyword_confidence.keywords()):
    ratio, p = keyword_confidence.keyword_confidence(key, condition)
    dd.append(
        {
            "keyword": key,
            "odds_ratio": ratio,
            "p_value": p,
        }
    )

# %%
v = pd.DataFrame(dd)
v.sort_values("p_value", ascending=True)

# %%
