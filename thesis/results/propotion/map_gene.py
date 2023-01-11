# %%
from src.plot.util.process_report import count_gene_type
import pandas as pd
import os
from dotenv import load_dotenv
from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics.condition import ConditionEq, ConditionAnd

load_dotenv()

THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]
TB_SAVEDIR = os.path.join(THESIS_TB_PATH, "results")


if not os.path.exists(TB_SAVEDIR):
    os.makedirs(TB_SAVEDIR)

# %%


def squash_columns(df: pd.DataFrame, threshold: int) -> pd.DataFrame:
    """値が少ないカラムはotherでまとめる"""
    df.fillna(0, inplace=True)
    columns_exclude = list(df.max().loc[lambda x: x <= threshold].index)  # type: ignore
    result = df.loc[:, ~df.columns.isin(columns_exclude)].copy()  # type: ignore
    print("squashed columns: ", columns_exclude)
    result["others"] = list(df[columns_exclude].sum(axis=1))
    return result


# %%

N = 10
for biosample in ["HepG2", "K562"]:
    report = SampleSetECLIP(ConditionEq("Biosample name", biosample)).report
    data = squash_columns(count_gene_type(report), 0)
    seri = data.sum().sort_values(ascending=False)
    seri.tail(seri.shape[0] - N + 1).sum()
    tmp = seri.head(N - 1).append(
        pd.Series({"others": seri.tail(seri.shape[0] - N + 1).sum()})
    )
    ratio_table = pd.DataFrame({"count": tmp, "proportion": tmp / tmp.sum() * 100})
    ratio_table.style.format(
        {"count": "{:,.0f}", "proportion": "{:.2f}\%"}, escape="latex"
    ).format_index(escape="latex", axis=0).format_index(
        escape="latex", axis=1
    ).to_latex(
        os.path.join(TB_SAVEDIR, f"gene_type_{biosample}.tex"),
        column_format="lrr",
        position="tbp",
        position_float="centering",
        hrules=True,
        caption=f"bedtools intersectによって変換された後の遺伝子の件数と割合 ({biosample})",
        label=f"tab:gene_type_{biosample}",
    )
# %%

# DKC1 HepG2, TROVE2 K562
targets = {"HepG2": "DKC1", "K562": "TROVE2"}
for biosample, target in targets.items():
    report = SampleSetECLIP(
        ConditionAnd(
            [
                ConditionEq("Biosample name", biosample),
                ConditionEq("Target label", target),
            ]
        )
    ).report
    data = squash_columns(count_gene_type(report), 0).iloc[0]
    ratio_table = pd.DataFrame(
        {"count": data, "proportion": data / data.sum() * 100}
    ).sort_values("count", ascending=False)
    print(biosample, target)
    print(ratio_table.head())

# %%
# 90%以上がprotein_codingのものを数える

for biosample in ["HepG2", "K562"]:
    report = SampleSetECLIP(ConditionEq("Biosample name", biosample)).report
    data = squash_columns(count_gene_type(report), 0)
    total_value = data.sum(axis=1).to_numpy().reshape(-1, 1)
    data = data.div(total_value, axis=1)
    result = (data["protein_coding"] >= 0.9).sum()
    print(data.shape, result)
# %%
