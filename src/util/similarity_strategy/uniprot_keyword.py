from src.util.similarity_strategy.interface import ProteinSimilarityStrategy
from src.eclip.uniprot.keyword import Keyword
from typing import List, Dict
import math
import pandas as pd
from tqdm import tqdm


class KeywordCosine(ProteinSimilarityStrategy):
    def _protein_similarity(self) -> pd.DataFrame:
        from scipy.spatial.distance import cosine, pdist, squareform
        from sklearn.feature_extraction.text import TfidfVectorizer

        from src.util.uniprot import load_keyword_report

        keywords = load_keyword_report(False)

        keywords["Keyword ID"] = keywords["Keyword ID"].map(
            lambda x: [key.strip() for key in x.split(";")]
        )
        keywords["label"] = "sp|" + keywords["Entry"] + "|" + keywords["Entry Name"]
        vectorizer = TfidfVectorizer()
        print(keywords["Keyword ID"].head())
        X = vectorizer.fit_transform(
            keywords["Keyword ID"].map(lambda x: " ".join([y[2:] for y in x]))
        )
        print(vectorizer.get_feature_names_out())
        return pd.DataFrame(
            squareform(pdist(X.toarray(), cosine)),  # type: ignore
            index=keywords["label"].tolist(),
            columns=keywords["label"].tolist(),
        )

    def __repr__(self) -> str:
        return "KeywordCosine"

    @property
    def lower_better(self) -> bool:
        return True


class KeywordAA(ProteinSimilarityStrategy):
    def _protein_similarity(self) -> pd.DataFrame:
        self.__load()
        proteins = self._edges["Protein"].unique()
        data = pd.DataFrame(columns=proteins, index=proteins, dtype=int)
        for protein_a in tqdm(proteins):
            for protein_b in proteins:
                # https://networkx.org/documentation/networkx-1.9/_modules/networkx/algorithms/link_prediction.html#adamic_adar_index
                data.loc[protein_a, protein_b] = sum(
                    1 / math.log(self.keyword_degree(keyword))
                    for keyword in self.common_keywords(protein_a, protein_b)
                )
        return data

    def __load(self):
        # 同じタンパク質でKeywordAAを計算するとゼロ除算が発生する
        # 1回しか出てこないKeywordを省く
        print("load keyword report")
        self._edges = Keyword.protein_keyword_table()
        keyword_counts = self._edges["Keyword"].value_counts()
        ignore_keywords = keyword_counts[keyword_counts == 1].index
        self._edges = self._edges[~self._edges["Keyword"].isin(ignore_keywords)]
        self._protein2keyword = self._edges.groupby("Protein")["Keyword"].apply(list)
        self._keyword2protein = self._edges.groupby("Keyword")["Protein"].apply(list)

    def common_keywords(self, protein_a: str, protein_b: str) -> List:
        if not hasattr(self, "_protein2keyword"):
            self.__load()
        return list(
            set(self._protein2keyword[protein_a])
            & set(self._protein2keyword[protein_b])
        )

    def keyword_degree(self, keyword) -> float:
        if not hasattr(self, "_keyword2protein"):
            self.__load()
        return len(self._keyword2protein[keyword])

    def __repr__(self) -> str:
        return "KeywordAA"

    @property
    def lower_better(self) -> bool:
        return False

    def _idmapping(self) -> Dict[str, str]:
        if not hasattr(self, "_edges"):
            self.__load()
        proteins = self._edges["Protein"].unique()
        return {protein: protein for protein in proteins}
