def test_tape_cosine():
    from src.util.tape import tape_cosine_similarity, load_embedding
    from scipy.spatial.distance import cosine
    import numpy as np

    df = tape_cosine_similarity()
    embedding = load_embedding()
    assert df.shape[0] == df.shape[1]
    assert df.shape[0] == len(embedding)
    for index, row in df.head(10).iterrows():
        for col, value in row.items():
            expected = cosine(embedding[str(index)], embedding[str(col)])
            assert np.isclose(value, expected, atol=1e-7)
