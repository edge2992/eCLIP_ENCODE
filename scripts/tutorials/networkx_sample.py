# %%
import networkx as nx
import matplotlib.pyplot as plt

# %%
G = nx.Graph()
G2 = nx.Graph([(0, 1), (1, 2), (2, 3), (3, 0)])

# %%
nx.draw(G2, with_labels=True)
plt.show()

# %%
