# keywordごとにstringdbのスコアの分布を見てみる

# %%
import os

import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

from src.plot.foldseek_stringdb_investigation.plot_utils import (
    construct_plotting_data,
    plot_boxplot_by_keyword,
)
from src.plot.interaction_metrics.representative import (
    convert_to_dict_exp_pair_by_keyword,
    target_report,
)
from src.util.metrics import Metrics
from src.util.similarity_strategy import TAPE, DirectStringScore, Jaccard, Simpson

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
keyword_experiment_pair = convert_to_dict_exp_pair_by_keyword(data)
# %%

investigation_metrics = [
    "stringdb_score",
    "stringdb_ascore",
    "stringdb_escore",
    "stringdb_tscore",
]
plot_data = construct_plotting_data(
    data, keyword_experiment_pair, investigation_metrics
)

for met in investigation_metrics:
    fig = plot_boxplot_by_keyword(plot_data, met)
    fig.savefig(os.path.join(SAVEDIR, f"boxplot_{met}.png"), bbox_inches="tight")

# %%
