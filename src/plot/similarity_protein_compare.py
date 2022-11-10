# タンパク質同士の類似度について
# 1. MSAの距離
# 2. TAPEのコサイン類似度
# 3. keywordのtfidfのコサイン類似度
# を比較する
# %%
import pandas as pd
import seaborn as sns
import os
from dotenv import load_dotenv

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

from src.util.similarity_protein import Similarity
from src.util.similarity_strategy import MSA, TAPE, KeywordCosine

# %%
similarity = Similarity()

similarity.setStrategy(MSA())
msa = similarity.executeStrategy()

similarity.setStrategy(TAPE())
tape = similarity.executeStrategy()

similarity.setStrategy(KeywordCosine())
keyword = similarity.executeStrategy()

assert msa.shape == tape.shape == keyword.shape
assert (msa.columns == tape.columns).all()
assert (keyword.columns == tape.columns).all()

# %%

data = pd.DataFrame(
    {
        "MSA": msa.to_numpy().reshape(-1),
        "TAPE": tape.to_numpy().reshape(-1),
        "Keyword": keyword.to_numpy().reshape(-1),
    }
)

splot = sns.pairplot(data, plot_kws={"alpha": 0.1})
splot.fig.suptitle("Similarity of proteins")
splot.fig.subplots_adjust(top=0.9)
splot.fig.savefig(
    os.path.join(PROJECT_PATH, "src/plot/img", "similarity_protein_compare.png")
)

# %%
