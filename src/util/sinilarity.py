def similarity_set(gene_list1, gene_list2):
    """二つの遺伝子の集合の類似度を計算する"""
    # TODO: 実装方針を立てる
    return len(set(gene_list1) & set(gene_list2))


def load_protein_sequence_similarity():
    # clustal Omegaのpercent identity Matrixを読み取って、pd.DataFrameにする
    pass
