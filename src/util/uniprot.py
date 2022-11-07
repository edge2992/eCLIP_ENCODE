import pandas as pd
import os
from dotenv import load_dotenv

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]


def idmapping():
    """uniprot idmappingからダウンロードしたテーブルデータを使用して、uniprotのproteinidとeCLIPで使われているタンパク質の対応表を作成する"""
    df = (
        pd.read_table(os.path.join(PROJECT_PATH, "data/uniprot", "reviewed.tsv"))
        .sort_values("Length", ascending=False)
        .drop_duplicates("From")
    )
    df["label"] = "sp|" + df["Entry"] + "|" + df["Entry Name"]
    return df[["From", "label"]].set_index("From").to_dict()["label"]
