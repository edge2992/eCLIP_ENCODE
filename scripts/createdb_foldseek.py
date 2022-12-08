# %%

from src.util.uniprot import idmapping
import os
import shutil
import pandas as pd

# %%

ALDB_DIR = "/home/edge2992/H/MYWORK/alphafolddb"
RDB_STRUCT_DIR = "/home/edge2992/H/MYWORK/eCLIP_ENCODE/data/afdb/rbp"


def idmapping_mmcif():
    uniprot_dict = idmapping()
    result = {}
    for key, value in uniprot_dict.items():
        result[key] = "AF-{}-F1-model_v3.cif.gz".format(value.split("|")[1])
    return result


# %%
# ファイルの存在確認をする
data = idmapping_mmcif()

if not os.path.exists(RDB_STRUCT_DIR):
    os.makedirs(RDB_STRUCT_DIR)

for key, value in data.items():
    target_path = os.path.join(ALDB_DIR, value)
    if not os.path.exists(target_path):
        print("[not exist] ", key, value)
        pass
    shutil.copy2(target_path, RDB_STRUCT_DIR + "/")


# %%
# foldseek createdb rdb queryDB
# foldseek search queryDB targetDB aln tmpFolder -a
# foldseek aln2tmscore queryDB targetDB aln aln_tmscore
# foldseek createtsv queryDB targetDB aln_tmscore aln_tmscore.tsv

# %%

FOLDSEEK_PATH = "/home/edge2992/H/MYWORK/eCLIP_ENCODE/data/afdb/aln_tmscore.tsv"
ALN_TMSCORE_COLUMNS = [
    "query",
    "target",
    "TMscore",
    "translation_0",
    "translation_1",
    "translation_2",
    "rotation_0_0",
    "rotation_0_1",
    "rotation_0_2",
    "rotation_1_0",
    "rotation_1_1",
    "rotation_1_2",
    "rotation_2_0",
    "rotation_2_1",
    "rotation_2_2",
]
df = pd.read_table(
    FOLDSEEK_PATH, names=ALN_TMSCORE_COLUMNS, sep="[ \t]", engine="python"
)
df.head()

protein_list = []
convert_dict = {v: k for k, v in idmapping_mmcif().items()}

for index, row in df.iterrows():
    protein_list.append(
        {
            "query_protein": convert_dict[row["query"]],
            "target_protein": convert_dict[row["target"]],
        }
    )

sample = pd.concat([pd.DataFrame(protein_list), df], axis=1)
sample.head()

# %%
