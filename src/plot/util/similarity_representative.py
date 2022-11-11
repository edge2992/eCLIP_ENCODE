from typing import List, Tuple
import pandas as pd
import numpy as np
import os
from dotenv import load_dotenv
from scipy.spatial.distance import squareform

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

from src.util.similarity_protein import ProteinSimilarity
from src.util.similarity_strategy import TAPE, KeywordCosine, MSA
from src.util.uniprot import load_keyword_report


def get_ranked_similarity() -> Tuple[pd.DataFrame, List[str]]:
    similarity = ProteinSimilarity()

    def get_flatten_tri(strategy):
        similarity.setStrategy(strategy)
        return similarity.flatten_tri(similarity.executeStrategy(), False)

    def get_protein_column(strategy):
        similarity.setStrategy(strategy)
        return list(similarity.executeStrategy().columns)

    data = pd.DataFrame(
        {
            "TAPE": get_flatten_tri(TAPE()),
            "Keyword": get_flatten_tri(KeywordCosine()),
            "MSA": get_flatten_tri(MSA()),
        }
    )
    data["rank"] = (data["TAPE"] + data["Keyword"]).rank()
    return data, get_protein_column(TAPE())


def get_topN_protein_pair(topN=100):
    data, protein_columns = get_ranked_similarity()
    mat = np.triu(squareform(data["rank"] <= topN))
    index_topN = np.where(mat == True)
    protein_data = pd.DataFrame(
        {
            "protein1": [protein_columns[i] for i in index_topN[0]],
            "protein2": [protein_columns[i] for i in index_topN[1]],
            "value": squareform(data["TAPE"] + data["Keyword"])[index_topN],
        }
    )
    return protein_data.sort_values("value").reset_index(drop=True)


def get_topN_protein_pair_keyword_appended(topN=100):
    return pd.merge(
        pd.merge(
            get_topN_protein_pair(topN),
            load_keyword_report()[["From", "Keywords"]]
            .rename({"Keywords": "Keywords1"}, axis=1)
            .drop_duplicates("From"),
            left_on="protein1",
            right_on="From",
        ),
        load_keyword_report()[["From", "Keywords"]]
        .rename({"Keywords": "Keywords2"}, axis=1)
        .drop_duplicates("From"),
        left_on="protein2",
        right_on="From",
    )[["protein1", "protein2", "Keywords1", "Keywords2", "value"]]
