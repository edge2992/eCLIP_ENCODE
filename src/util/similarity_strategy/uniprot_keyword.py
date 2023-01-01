from src.util.similarity_strategy.interface import ProteinSimilarityStrategy
import pandas as pd


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
