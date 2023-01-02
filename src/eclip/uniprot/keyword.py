from typing import List
import pandas as pd

from src.util.uniprot import load_keyword_report


class Singleton(object):
    def __new__(cls, *args, **kargs):
        if not hasattr(cls, "_instance"):
            cls._instance = super(Singleton, cls).__new__(cls)
        return cls._instance


class Keyword(Singleton):
    def __init__(self):
        if not hasattr(self, "_keywords"):
            self._keywords = self.__load()

    def __call__(self, protein) -> List[str]:
        keywords = self._keywords.loc[protein, "Keywords"]
        assert isinstance(keywords, str)
        return [key.strip() for key in keywords.split(";")]

    @classmethod
    def __load(cls):
        print("load keyword report")
        return load_keyword_report().set_index("From")

    @classmethod
    def keywords(cls) -> List[str]:
        """keywordを全て取得する"""
        keyword_list = []
        for _, row in cls()._keywords.iterrows():
            keyword_list.extend([key.strip() for key in row["Keywords"].split(";")])
        return list(set(keyword_list))

    @classmethod
    def protein_keyword_table(cls) -> pd.DataFrame:
        """proteinとkeywordの関係を表すDataFrameを返す"""
        return (
            cls()
            ._keywords.apply(lambda row: row["Keywords"].split(";"), axis=1)
            .reset_index()
            .rename(columns={"From": "Protein", 0: "Keyword"})
            .explode("Keyword")
            .reset_index(drop=True)
        )
