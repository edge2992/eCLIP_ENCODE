# %%
from src.eclip.uniprot.keyword import Keyword

keyword = Keyword()
# %%
len(keyword.keywords())
# %%
len(keyword.keywords())


# %%
df = keyword.protein_keyword_table()

# %%

df.nunique()
# %%
df["Protein"].value_counts()

# %%
df["Keyword"].value_counts()

# %%
