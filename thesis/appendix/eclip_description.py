# %%
import os

from dotenv import load_dotenv

from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics.condition import ConditionEq

load_dotenv()

THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]

SAVEDIR = os.path.join(THESIS_TB_PATH, "appendix")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# %%

reports = {
    k: SampleSetECLIP(ConditionEq("Biosample name", k)).report
    for k in ["HepG2", "K562", "adrenal gland"]
}

# %%
for biosample in ["HepG2", "K562", "adrenal gland"]:
    data = reports[biosample][["Dataset", "Accession", "Target label"]].copy()
    data["Dataset"] = data["Dataset"].apply(lambda x: x.split("/")[2])
    print(biosample, data.shape)
    data.rename({"Accession": "peaks"}, axis=1).to_latex(
        os.path.join(SAVEDIR, f"encode_{biosample}.tex"),
        column_format="llll",
        longtable=True,
        caption=f"Biosampleが{biosample}であるデータセットと解析ファイルのAccession番号",
        label=f"tab:encode_{biosample}",
        index=False,
    )

# %%
