# %%
import pandas as pd


COMPREHENSIVE_PATH = "~/Resource/gencode.v24.annotation.gtf"
GENE_ANOTATION_PATH = "~/Resource/gencode.v24.basic.annotation.gtf"
LNCRNA_ANOTATION_PATH = "~/Resource/gencode.v24.long_noncoding_RNAs.gtf"

# %%

df = pd.read_table(COMPREHENSIVE_PATH, skiprows=5, header=None)

# %%
df.head()
print(df[2].value_counts())

# %%
gene = pd.read_table(GENE_ANOTATION_PATH, skiprows=5, header=None)
lncrna = pd.read_table(LNCRNA_ANOTATION_PATH, skiprows=5, header=None)

# %%
print(gene[2].value_counts())

# %%
print(lncrna[2].value_counts())

# %%
