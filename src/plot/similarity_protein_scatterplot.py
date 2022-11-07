# %%
# 配列類似度とlift値の散布図を描画する
from typing import List
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from dotenv import load_dotenv
from src.util.similarity import lift_protein, load_protein_sequence_similarity
from src.util.bedfile import load_replicateIDR_report

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]


# %%
def format_similarity(column: List[str]):
    sequence_similarity = load_protein_sequence_similarity()
    data = pd.DataFrame(index=column, columns=column)
    for index, row in data.iterrows():
        for col, _ in row.items():
            data.loc[index, col] = sequence_similarity.loc[  # type: ignore
                str(index).split()[0], str(col).split()[0]
            ]
    return data


lift = lift_protein(load_replicateIDR_report())
sequence_similarity = format_similarity(list(lift.columns))

# %%
assert sequence_similarity.shape == lift.shape
assert (sequence_similarity.columns == lift.columns).all()
a = sequence_similarity.to_numpy().reshape(-1)
b = lift.to_numpy().reshape(-1)
assert len(a) == len(b)

# %%

fig, ax = plt.subplots(figsize=(10, 10))
sns.scatterplot(x=a, y=b, ax=ax)
ax.set_yscale("log")
ax.set_xlabel("<- sequence distance calculated from Clustal Omega", fontsize=20)
ax.set_ylabel("lift ->", fontsize=20)
fig.savefig(
    os.path.join(PROJECT_PATH, "src/plot/img", "similarity_protein_scatter.png")
)

# %%
