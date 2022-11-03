import numpy as np
import pandas as pd
from typing import List


def split_list(lst: List[str], n: int) -> List[List[str]]:
    """Split the list into N pieces by splitting it from the front"""
    return [list(tmp) for tmp in np.array_split(lst, n)]


def split_dataframe(df: pd.DataFrame, n: int) -> List[pd.DataFrame]:
    """Split the dataframe into N pieces by splitting it from the front"""
    return [df.loc[tmp] for tmp in np.array_split(df.index, n)]
