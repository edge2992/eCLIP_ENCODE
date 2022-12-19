# %%
import os

import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv

from src.eclip.encodecondition import ECLIP_SAMPLESETS
from src.eclip.sampleset import SampleSetECLIP
from src.util.metrics import Metrics
from src.util.similarity_strategy import Jaccard, PeakStrategy

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
SAVEDIR = os.path.join(PROJECT_PATH, "src/plot/img", "peak_gene_jaccard")

if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


# %%
STRATEGIES = [
    Jaccard(),
    PeakStrategy(
        metrics="jaccard",
    ),
    PeakStrategy(metrics="union-intersection"),
    PeakStrategy(metrics="intersection"),
    PeakStrategy(metrics="n_intersections"),
]

for ss in ECLIP_SAMPLESETS:
    sampleset = SampleSetECLIP(ss)
    data = Metrics(sampleset.report)(STRATEGIES)
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.scatterplot(data=data, x="Jaccard", y="peak jaccard", ax=ax)
    ax.set_yscale("log")
    ax.set_title(str(sampleset))
    fig.savefig(os.path.join(SAVEDIR, f"{sampleset}.png"), bbox_inches="tight")

# %%
