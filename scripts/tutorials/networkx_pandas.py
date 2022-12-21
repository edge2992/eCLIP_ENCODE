# %%
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

# %%

r = np.random.RandomState(seed=5)
ints = r.random_integers(1, 10, size=(3, 2))
a = ["a", "b", "c"]
b = ["d", "a", "e"]
df = pd.DataFrame(ints, columns=["weight", "cost"])
df.head()
df[0] = a
df["b"] = b
df.head()
G = nx.from_pandas_edgelist(df, 0, "b", ["weight", "cost"])  # type: ignore
nx.draw(G, with_labels=True)
plt.plot()

# %%
# Edges
# https://www.roelpeters.be/python-networkx-set-node-attributes-from-pandas-dataframe/
df_edges = pd.DataFrame(
    {
        "source": ["John", "John", "Jane", "Alice", "Tom", "Helen", "William"],
        "target": ["Jane", "Helen", "Tom", "Tom", "Helen", "William", "Alice"],
        "years": [2, 3, 5, 1, 2, 8, 3],
    }
)
# Nodes
df_nodes = pd.DataFrame(
    {
        "Name": ["John", "Jane", "Alice", "Tom", "Helen", "William"],
        "Gender": ["Male", "Female", "Female", "Male", "Female", "Male"],
    }
)
node_colors = {"Male": "blue", "Female": "red"}
df_nodes["node_color"] = df_nodes["Gender"].map(node_colors)

G = nx.from_pandas_edgelist(
    df_edges, source="source", target="target", edge_attr="years"
)

node_attr = df_nodes.set_index("Name").to_dict(orient="index")
print(node_attr)
nx.set_node_attributes(G, node_attr)

plt.figure(figsize=(6, 6))
nx.draw_networkx(
    G,
    pos=nx.kamada_kawai_layout(G, weight="years"),
    node_color=[G.nodes[n]["node_color"] for n in G.nodes],
    width=[G.edges[e]["years"] for e in G.edges],
    with_labels=True,
)
plt.plot()

# %%
