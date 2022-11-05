from dotenv import load_dotenv
import os

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]


def test_uniprot_idmapping():
    import pandas as pd

    PROTEIN_TSV = os.path.join(PROJECT_PATH, "data", "uniprot", "reviewed.tsv")
    if not os.path.exists(PROTEIN_TSV):
        assert False, "file does not exist"
    df = pd.read_table(PROTEIN_TSV)

    # DDX21, TAF15がreviewedで二つある（それ以外は1つ)
    values = df["From"].value_counts()
    duplicated = list(values[values > 1].index)
    print(df[df["From"].isin(duplicated)])

    assert len(df["From"].unique()) == 159

    # for _, row in df.iterrows():
    #     assert row["From"] == row["Entry Name"].split("_")[0]
