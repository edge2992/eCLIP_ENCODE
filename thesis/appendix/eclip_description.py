# %%
from dotenv import load_dotenv
import os

from src.util.bedfile import load_replicateIDR_report

load_dotenv()

THESIS_TB_PATH = os.environ["THESIS_TB_PATH"]

SAVEDIR = os.path.join(THESIS_TB_PATH, "appendix")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# %%

df = load_replicateIDR_report()

# %%

df.head()

# %%
df.columns
# %%
encode_tb = df[["Dataset", "Accession", "Biosample name", "Target label"]].copy()
encode_tb["Dataset"] = encode_tb["Dataset"].apply(lambda x: x.split("/")[2])

# %%
encode_tb.head()

# %%
for label, grouped in encode_tb.groupby("Biosample name"):
    grouped.sort_values("Target label")[
        ["Dataset", "Accession", "Target label"]
    ].rename({"Accession": "peaks"}, axis=1).to_latex(
        os.path.join(SAVEDIR, f"encode_{label}.tex"),
        column_format="llll",
        longtable=True,
        caption=f"Biosampleが{label}であるデータセットと解析ファイルのAccession番号",
        label=f"tab:encode_{label}",
        index=False,
    )

# %%
