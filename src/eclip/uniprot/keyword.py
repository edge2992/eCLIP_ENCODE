from typing import List

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
