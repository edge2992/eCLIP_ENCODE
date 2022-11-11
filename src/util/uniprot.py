import pandas as pd
import os
from dotenv import load_dotenv

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]
KEYWORD_FILE = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/keyword_list.tsv"
KEYWORD_PROTEIN = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_keyword.tsv"


def idmapping():
    """uniprot idmappingからダウンロードしたテーブルデータを使用して、uniprotのproteinidとeCLIPで使われているタンパク質の対応表を作成する"""
    df = (
        pd.read_table(os.path.join(PROJECT_PATH, "data/uniprot", "reviewed.tsv"))
        .sort_values("Length", ascending=False)
        .drop_duplicates("From")
    )
    df["label"] = "sp|" + df["Entry"] + "|" + df["Entry Name"]
    return df[["From", "label"]].set_index("From").to_dict()["label"]


def load_keyword_report(is_unique: bool = True) -> pd.DataFrame:
    """uniprotのkeywordを読み込む
    pd.DataFrame[From , Entry, Reviewed, Entry Name, Keywords, Keyword ID]

    is_unique: bool, default True trueの場合は、重複削除する
    """
    data = pd.read_table(KEYWORD_PROTEIN)
    if is_unique:
        return data.drop_duplicates("From", keep="first")
    return data


def load_keyword_list():
    """uniprotで使われているkeywordの一覧を返す"""
    return pd.read_table(KEYWORD_FILE).drop("Gene Ontologies", axis=1)


def load_splited_keyword_report(is_unique: bool = True) -> pd.DataFrame:
    """keywordを分割して読み込む"""
    tmp = (
        load_keyword_report(is_unique)
        .set_index("From")["Keyword ID"]
        .map(lambda x: [key.strip() for key in x.split(";")])
        .explode()
    )
    return pd.merge(tmp.reset_index(), load_keyword_list(), on="Keyword ID")


def keyword_count():
    """uniprotでkeywordの出現順と個数をまとめたDataFrameを返す"""
    counted: pd.Series = (
        pd.read_table(KEYWORD_PROTEIN)["Keyword ID"]
        .map(lambda x: [key.strip() for key in x.split(";")])
        .explode()
        .value_counts()
    ).rename("count")
    result = pd.concat(
        [counted, load_keyword_list().set_index("Keyword ID", drop=True)],
        axis=1,
        join="inner",
    )
    assert counted.shape[0] == result.shape[0]
    return result
