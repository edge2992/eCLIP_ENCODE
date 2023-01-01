from src.util.similarity_strategy import KeywordAA
from src.util.similarity_protein import ProteinSimilarity


def test_keywordAA_protein():
    """タンパク質のKeywordAAテーブルが作成できる"""
    import math

    handler = ProteinSimilarity()
    handler.setStrategy(KeywordAA())
    data = handler.executeStrategy()

    strategy = KeywordAA()
    for index, row in data.head().iterrows():
        for column, value in row.items():
            expcted = sum(
                1 / math.log(strategy.keyword_degree(keyword))
                for keyword in strategy.common_keywords(str(index), str(column))
            )
            assert value == expcted


def test_common_keywords():
    strategy = KeywordAA()
    protein_a = "IGF2BP3"
    protein_b = "YWHAG"
    common_key_list = strategy.common_keywords(protein_a, protein_b)
    assert len(common_key_list) == 4
    for keyword in common_key_list:
        print(keyword, strategy.keyword_degree(keyword))


def test_keywordAA_interaction():
    """KeywordAAテーブルを実験データごとのテーブルに変換できる"""
    pass
