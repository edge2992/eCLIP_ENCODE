import pandas as pd
import os
from dotenv import load_dotenv

load_dotenv()
GENE_GTF_PATH = os.environ["GENE_GTF_PATH"]


COLUMN_ENSEMBL_GTF = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]


def load_gtf(filename: str = GENE_GTF_PATH) -> pd.DataFrame:
    return pd.read_table(filename, names=COLUMN_ENSEMBL_GTF)


def gene_ids_gtf(filename: str = GENE_GTF_PATH):
    return (
        load_gtf(filename)["attribute"]
        .str.split('gene_id "', expand=True)[1]
        .str.split('";', expand=True)[0]
        .unique()
        .tolist()
    )
