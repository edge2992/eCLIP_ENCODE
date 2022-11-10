# RBPのRNAとの相互作用の類似度について
# 1.

# %%
import pandas as pd
import seaborn as sns
import os
from dotenv import load_dotenv

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

from src.util.similarity_protein import Similarity
from src.util.similarity_strategy import Lift, Dice, Jaccard, Simpson, Cosine

# %%
similarity = Similarity()
similarity.setStrategy(Lift())
lift = similarity.executeStrategy()

similarity.setStrategy(Dice())
dice = similarity.executeStrategy()

similarity.setStrategy(Jaccard())
jaccard = similarity.executeStrategy()

similarity.setStrategy(Simpson())
simpson = similarity.executeStrategy()

similarity.setStrategy(Cosine())
cosine = similarity.executeStrategy()

# %%

data = pd.DataFrame(
    {
        "Lift": lift.to_numpy().reshape(-1),
        "Dice": dice.to_numpy().reshape(-1),
        "Jaccard": jaccard.to_numpy().reshape(-1),
        "Simpson": simpson.to_numpy().reshape(-1),
        "Cosine": cosine.to_numpy().reshape(-1),
    }
)
splot = sns.pairplot(data, plot_kws={"alpha": 0.1})
splot.fig.suptitle("Similarity of proteins")
splot.fig.subplots_adjust(top=0.9)
splot.fig.savefig(
    os.path.join(PROJECT_PATH, "src/plot/img", "similarity_interaction_compare.png")
)

# %%

# TODO: タンパク質の類似度指標で色付けする
# 1. 連続値で色付け
# 2. TAPEとKeywordのsimilarityの上位10%を色付け
# 2.1 上位10%がなにかを出力で示す
