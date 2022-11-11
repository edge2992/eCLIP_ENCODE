# RBPのRNAとの相互作用の類似度について
# 1.

# %%
import os

import pandas as pd
import seaborn as sns
from dotenv import load_dotenv

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

from src.util.similarity_protein import InteractionSimilarity, ProteinSimilarity
from src.util.similarity_strategy import (
    MSA,
    TAPE,
    Cosine,
    Dice,
    Jaccard,
    KeywordCosine,
    Lift,
    Simpson,
)

# %%
similarity = InteractionSimilarity()
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
        "Lift": similarity.flatten_tri(lift),
        "Dice": similarity.flatten_tri(dice),
        "Jaccard": similarity.flatten_tri(jaccard),
        "Simpson": similarity.flatten_tri(simpson),
        "Cosine": similarity.flatten_tri(cosine),
    }
)
splot = sns.pairplot(data, plot_kws={"alpha": 0.1})
splot.fig.suptitle("Similarity of proteins")
splot.fig.subplots_adjust(top=0.9)
for ax in splot.axes.flatten():
    if ax.get_xlabel() == "Lift":
        ax.set(xscale="log")
    if ax.get_ylabel() == "Lift":
        ax.set(yscale="log")
splot.fig.savefig(
    os.path.join(
        PROJECT_PATH,
        "src/plot/img",
        "similarity_interaction_compare",
        "interaction.png",
    )
)

# %%
similarity = ProteinSimilarity()
similarity.setStrategy(TAPE())
tape = similarity.executeStrategy()
tape = similarity.strategy.transform(tape)
similarity.setStrategy(MSA())
msa = similarity.executeStrategy()
msa = similarity.strategy.transform(msa)
similarity.setStrategy(KeywordCosine())
keyword = similarity.executeStrategy()
keyword = similarity.strategy.transform(keyword)

# %%
data_protein_interaction = pd.DataFrame(
    {
        "Lift": similarity.flatten_tri(lift, False),
        "Dice": similarity.flatten_tri(dice, False),
        "Jaccard": similarity.flatten_tri(jaccard, False),
        "Simpson": similarity.flatten_tri(simpson, False),
        "Cosine": similarity.flatten_tri(cosine, False),
        "TAPE": similarity.flatten_tri(tape, False),
        "MSA": similarity.flatten_tri(msa, False),
        "Keyword": similarity.flatten_tri(keyword, False),
    }
)

splot = sns.pairplot(
    data_protein_interaction.sample(3000),
    x_vars=["Lift", "Dice", "Jaccard", "Simpson", "Cosine"],
    y_vars=["TAPE", "MSA", "Keyword"],
    plot_kws={"alpha": 0.1},
)
for ax in splot.axes.flatten():
    if ax.get_xlabel() == "Lift":
        ax.set(xscale="log")
splot.fig.savefig(
    os.path.join(
        PROJECT_PATH,
        "src/plot/img",
        "similarity_interaction_compare",
        "protein_interaction.png",
    )
)

# %%
# タンパク質の上位N件のmatrixをDatasetのmatrixに変換する
from src.plot.util.similarity_representative import get_ranked_similarity
from scipy.spatial.distance import squareform

data, protein_column = get_ranked_similarity()

# %%
TOPN = 100
(data["rank"] <= TOPN).to_numpy()
similarity = ProteinSimilarity()
similarity.setStrategy(TAPE())
representative = similarity.strategy.transform(
    pd.DataFrame(
        squareform(data["rank"] <= TOPN), columns=protein_column, index=protein_column
    )
)
# %%
data_protein_interaction["representative"] = similarity.flatten_tri(
    representative, False
)
data_protein_interaction["representative"] = data_protein_interaction[
    "representative"
].astype(bool)

plot_data = pd.concat(
    [
        data_protein_interaction[data_protein_interaction["representative"]],
        data_protein_interaction[~data_protein_interaction["representative"]].sample(
            3000
        ),
    ]
)


# %%
data_protein_interaction.head()
splot = sns.pairplot(
    plot_data,
    x_vars=["Lift", "Dice", "Jaccard", "Simpson", "Cosine"],
    y_vars=["TAPE", "MSA", "Keyword"],
    hue="representative",
    plot_kws={"alpha": 0.1},
)
for ax in splot.axes.flatten():
    if ax.get_xlabel() == "Lift":
        ax.set(xscale="log")
splot.fig.suptitle(
    "Similarity of proteins. Top 100 protein pairs + 3000 random pairs. (from 12561 pairs)"
)
splot.fig.supxlabel("Interaction similarity")
splot.fig.supylabel("Protein similarity")
splot.fig.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.9)
splot.fig.savefig(
    os.path.join(
        PROJECT_PATH,
        "src/plot/img",
        "similarity_interaction_compare",
        "protein_interaction_top100.png",
    )
)
# %%
# distplotで比較する
# sns.displot(data_protein_interaction, x="Lift", hue="representative",  stat="density")

import matplotlib.pyplot as plt

interaction_similarity_label = ["Lift", "Dice", "Jaccard", "Simpson", "Cosine"]
fig, axes = plt.subplots(5, 1, figsize=(10, 15))
for ax, label in zip(axes, interaction_similarity_label):
    sns.kdeplot(
        data_protein_interaction,
        x=label,
        hue="representative",
        ax=ax,
        common_norm=False,
    )
    if label != "Lift":
        ax.set_xlim(-0.2, 1)
fig.subplots_adjust(hspace=0.4)
fig.savefig(
    os.path.join(
        PROJECT_PATH,
        "src/plot/img",
        "similarity_interaction_compare",
        "protein_interaction_top100_density.png",
    )
)
# %%
