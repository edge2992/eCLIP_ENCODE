# %%
import numpy as np
import matplotlib.pyplot as plt

plt.ion()

x = np.arange(-5, 5 + 1)
y = x**2


def graph(name):
    """2次関数をグラフ化
    Parameters
    ----------
    name : str
        作成したグラフの保存名
    """
    fig, ax = plt.subplots(1, 1)
    ax.set_title("title")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-5, 30)
    ax.plot(x, y, "o-", ms=4, label=r"$\alpha$")
    ax.plot(x, y + 2, "o-", ms=4, label=r"$\beta$")
    ax.plot(x, y + 3, "o-", ms=4, label=r"$\gamma$")
    ax.plot(x, y + 4, "o-", ms=4, label=r"$\dagger$")
    ax.legend()
    ax.grid(which="major", ls="-")
    fig.savefig(f"plotly-base_plt_{name}.png")


# この時点ではまだデフォルト設定
graph(name="before")

change_dct = {
    "font.size": 20,
    # "font.family": "Times New Roman",
    "mathtext.fontset": "stix",
    "figure.figsize": [16, 8],
    "figure.facecolor": (0, 0, 0, 0),
    "axes.facecolor": (0, 0, 0, 0),
    "savefig.facecolor": (0, 0, 0, 0),
    "figure.subplot.left": 0.11,
    "figure.subplot.right": 0.90,
    "figure.subplot.top": 0.95,
    "figure.subplot.bottom": 0.08,
    "figure.subplot.hspace": 0.1,
    "savefig.format": "pdf",
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "legend.edgecolor": "black",
    "legend.facecolor": "white",
    "legend.framealpha": 1,
}
# plt .rcParamsを変更
for key, val in change_dct.items():
    plt.rcParams[key] = val
# この時点ではデフォルトの値が変更されてる
graph(name="after")

# %%

import matplotlib.font_manager as fm
import pprint

font_list = [f.name for f in fm.fontManager.ttflist]

pprint.pprint(font_list)

# %%
import matplotlib

matplotlib.matplotlib_fname()

# %%
