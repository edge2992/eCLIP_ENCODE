def test_idmapping():
    from src.util.uniprot import idmapping
    from src.util.bedfile import load_replicateIDR_report

    eclip = load_replicateIDR_report()["Target label"].unique().tolist()
    d = idmapping()
    for protein in eclip:
        assert protein in d.keys()


def test_keyword_count():
    """keywordの出現順と個数をまとめたDataFrameを返す"""
    from src.util.uniprot import keyword_count

    df = keyword_count()
    print(df.head())
    print(df.columns)
    print(df.shape)


def test_load_keyword_reprot():
    """タンパク質ごとにキーワードをまとめたDataFrameを返す"""
    from src.util.uniprot import load_keyword_report

    df = load_keyword_report()
    print(df.head())
    print(df.columns)
    print(df.shape)


def test_load_splited_keyword_report():
    from src.util.uniprot import load_splited_keyword_report

    seri = load_splited_keyword_report()
    # assert seri.index.nunique() == 159
    assert seri["From"].nunique() == 159
    print(seri)
    print(seri.shape)
