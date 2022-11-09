from typing import Dict
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.spatial.distance import cosine

EMBEDDING_PATH = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed_tape.npz"


def load_embedding(how: str = "avg", filepath=EMBEDDING_PATH) -> Dict[str, np.ndarray]:
    """TAPE proteinのembeddingを読み込む

    Args:
        how (str, optional): avg or pool. Defaults to "avg".
        filepath (_type_, optional): ファイルパス. Defaults to EMBEDDING_PATH.

    Returns:
        Dict[str, np.ndarray]: 遺伝子名をキーとしたembeddingの辞書
    """

    arrays = np.load(filepath, allow_pickle=True)
    return {key: value[()][how] for key, value in arrays.items()}


def tape_cosine_similarity(how: str = "avg", filepath=EMBEDDING_PATH):
    """TAPEのcos類似度を計算する

    Args:
        how (str, optional): avg or pool. Defaults to "avg".
        filepath (_type_, optional): ファイルパス. Defaults to EMBEDDING_PATH.

    Returns:
        pd.DataFrame: cos類似度のデータフレーム
    """
    embedding = load_embedding(how=how, filepath=filepath)

    return pd.DataFrame(
        squareform(pdist(np.array(list(embedding.values())), cosine)),
        index=list(embedding.keys()),
        columns=list(embedding.keys()),
    )
