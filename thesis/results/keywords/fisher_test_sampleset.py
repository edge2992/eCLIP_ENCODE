# %%
import os
from tqdm import tqdm
from dotenv import load_dotenv
from src.util.metrics.keyword import KeywordConfidence
from src.util.similarity_strategy import Jaccard
from src.util.metrics import Metrics, ConditionGt, Condition
from thesis.utils.reportset import COMPARE_REPORT_SET
import pandas as pd

load_dotenv()

# THESIS_FIG_PATH = os.environ["THESIS_FIG_PATH"]
THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
# SAVEDIR = os.path.join(
#     THESIS_FIG_PATH, "results", "protein_interaction", "corr_heatmap"
# )
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results", "keywords")

if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)


# %%
def fisher_test_table(handler: KeywordConfidence, condition: Condition):
    table = []
    for keyword in tqdm(handler.keywords()):
        table.append(
            [
                keyword,
                *handler.keyword_confidence(keyword, condition),
            ]
        )
    return pd.DataFrame(
        table,
        columns=[
            "Keyword",
            "odds ratio",
            "p-value",
        ],
    )


# %%
for biosample in COMPARE_REPORT_SET:
    for key, sample_set in COMPARE_REPORT_SET[biosample].items():
        data = Metrics(sample_set.report)(Jaccard(), add_description=True)
        assert isinstance(data, pd.DataFrame)
        confidence_handler = KeywordConfidence(data)
        condition = ConditionGt(str(Jaccard()), 0.2)
        results = fisher_test_table(confidence_handler, condition)
        results.to_csv(
            os.path.join(TB_SAVEDIR, f"fisher_test_jaccard_{biosample}_{key}.csv"),
            index=False,
        )

# %%

for biosample in COMPARE_REPORT_SET:
    for key, sample_set in COMPARE_REPORT_SET[biosample].items():
        print(biosample, key)
        latex_key = "\_".join(key.split("_"))

        pd.read_csv(
            os.path.join(TB_SAVEDIR, f"fisher_test_jaccard_{biosample}_{key}.csv")
        ).sort_values("p-value").set_index("Keyword", drop=True).head(10).style.format(
            {"odds ratio": "{:.2f}", "p-value": "{:.2e}"}
        ).format_index(
            escape="latex", axis=0
        ).format_index(
            escape="latex", axis=0
        ).to_latex(
            os.path.join(TB_SAVEDIR, f"fisher_test_jaccard_{biosample}_{key}.tex"),
            column_format="lrr",
            position="htbp",
            position_float="centering",
            hrules=True,
            caption=f"キーワードに対するFisherの正確性検定の結果 ({latex_key} ({biosample})). p-valueの上位10件を表示した。",
            label=f"tab:jaccard_fisher_test_{biosample}_{key}",
        )


# %%
# %%

for biosample in COMPARE_REPORT_SET:
    for key, sample_set in COMPARE_REPORT_SET[biosample].items():
        t = pd.read_csv(
            os.path.join(TB_SAVEDIR, f"fisher_test_jaccard_{biosample}_{key}.csv")
        )
        print(t.head())
        print(t.shape)

# %%
