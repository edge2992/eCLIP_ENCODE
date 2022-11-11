import numpy as np


def test_past_load_label():
    """report.txtとsimilarity_matrixのラベル (タンパク質) が一致する"""
    """load_protein_sequence_similarityのテスト"""
    from src.util.similarity import load_protein_sequence_similarity
    from src.util.bedfile import load_replicateIDR_report

    df = load_protein_sequence_similarity()
    expected = load_replicateIDR_report()["Target label"].unique().tolist()
    assert df.shape[0] == df.shape[1]
    assert len(df.columns) == len(expected)
    for protein in expected:
        assert protein in df.columns
        assert protein in df.index
    data = df.to_numpy()
    assert (data == data.T).reshape(-1).all(), "対称行列である"


def test_past_lift_protein():
    """タンパク質のリフト値を計算する"""
    from src.util.similarity import lift_protein
    from src.util.bedfile import load_replicateIDR_report
    from src.plot.util.process_report import label_protein_biosample, gene_ids_eCLIP
    from src.util.bedfile import read_annotated_bed
    from src.util.get_bed_path import get_formatted_file_path
    from src.util.bed_format_strategy import FormatStrategy
    import pandas as pd

    def get_gene_list(report: pd.DataFrame, label_pb: str):
        seri = report[
            (report["Target label"] == label_pb.split()[0])
            & (report["Biosample name"] == label_pb.split()[1])
        ]
        assert len(seri) == 1
        assert label_protein_biosample(seri).values[0] == label_pb
        return (
            read_annotated_bed(
                get_formatted_file_path(seri.iloc[0], FormatStrategy.MAX)
            )["gene_id"]
            .dropna()
            .unique()
            .tolist()
        )

    N_TEST = 20
    report = load_replicateIDR_report().head(N_TEST)
    df = lift_protein(report)
    all_gene_count = len(gene_ids_eCLIP(report))
    assert df.shape[0] == df.shape[1]
    for index, row in df.iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(report, str(index))
            gene2 = get_gene_list(report, str(col))
            expected = (
                len(set(gene1) & set(gene2))
                / (len(gene1) * len(gene2))
                * all_gene_count
            )
            assert abs(value - expected) < 1e-9, f"{index}, {col} invalid"


def test_load_label():
    """report.txtとsimilarity_matrixのラベル (タンパク質) が一致する"""
    """Multiple Sequence Analysis Distanceのテスト"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.similarity_strategy import MSA

    similarity = ProteinSimilarity()
    similarity.setStrategy(MSA())

    df = similarity.executeStrategy()

    expected = load_replicateIDR_report()["Target label"].unique().tolist()
    assert df.shape[0] == df.shape[1]
    assert len(df.columns) == len(expected)
    for protein in expected:
        assert protein in df.columns
        assert protein in df.index
    data = df.to_numpy()
    assert (data == data.T).reshape(-1).all(), "対称行列である"


def test_tape():
    """report.txtとsimilarity_matrixのラベル (タンパク質) が一致する"""
    """Multiple Sequence Analysis Distanceのテスト"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.similarity_strategy import TAPE

    similarity = ProteinSimilarity()
    similarity.setStrategy(TAPE())

    df = similarity.executeStrategy()

    expected = load_replicateIDR_report()["Target label"].unique().tolist()
    assert df.shape[0] == df.shape[1]
    assert len(df.columns) == len(expected)
    for protein in expected:
        assert protein in df.columns
        assert protein in df.index
    data = df.to_numpy()
    assert (data == data.T).reshape(-1).all(), "対称行列である"


def test_lift_protein():
    """タンパク質のリフト値を計算する"""
    from src.util.bedfile import load_replicateIDR_report
    from src.plot.util.process_report import gene_ids_eCLIP
    from src.util.similarity_protein import InteractionSimilarity
    from src.util.similarity_strategy import Lift
    from tests.util.similarity import get_gene_list

    N_TEST = 20
    report = load_replicateIDR_report().head(N_TEST)
    all_gene_count = len(gene_ids_eCLIP(report))

    similarity = InteractionSimilarity()
    similarity.setStrategy(Lift(report=report, label_method=lambda df: df["Accession"]))
    df = similarity.executeStrategy()
    assert df.shape[0] == df.shape[1]
    for index, row in df.iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(report, str(index))
            gene2 = get_gene_list(report, str(col))
            expected = (
                len(set(gene1) & set(gene2))
                / (len(gene1) * len(gene2))
                * all_gene_count
            )
            assert np.isclose(value, expected, atol=1e-9), f"{index}, {col} invalid"


def test_dice_protein():
    """タンパク質のdice値を計算する"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import InteractionSimilarity
    from src.util.similarity_strategy import Dice
    from tests.util.similarity import get_gene_list

    N_TEST = 20
    report = load_replicateIDR_report().head(N_TEST)

    similarity = InteractionSimilarity()
    similarity.setStrategy(Dice(report=report, label_method=lambda df: df["Accession"]))
    df = similarity.executeStrategy()
    assert df.shape[0] == df.shape[1]
    for index, row in df.iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(report, str(index))
            gene2 = get_gene_list(report, str(col))
            expected = len(set(gene1) & set(gene2)) * 2 / (len(gene1) + len(gene2))
            assert np.isclose(value, expected, atol=1e-9), f"{index}, {col} invalid"


def test_jaccard_protein():
    """タンパク質のjaccard値を計算する"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import InteractionSimilarity
    from src.util.similarity_strategy import Jaccard
    from tests.util.similarity import get_gene_list

    N_TEST = 20
    report = load_replicateIDR_report().head(N_TEST)

    similarity = InteractionSimilarity()
    similarity.setStrategy(
        Jaccard(report=report, label_method=lambda df: df["Accession"])
    )
    df = similarity.executeStrategy()
    assert df.shape[0] == df.shape[1]
    for index, row in df.iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(report, str(index))
            gene2 = get_gene_list(report, str(col))
            expected = len(set(gene1) & set(gene2)) / len(set(gene1) | set(gene2))
            assert np.isclose(value, expected, atol=1e-9), f"{index}, {col} invalid"


def test_simpson_protein():
    """タンパク質のsimpson値を計算する"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import InteractionSimilarity
    from src.util.similarity_strategy import Simpson
    from tests.util.similarity import get_gene_list

    N_TEST = 20
    report = load_replicateIDR_report().head(N_TEST)

    similarity = InteractionSimilarity()
    similarity.setStrategy(
        Simpson(report=report, label_method=lambda df: df["Accession"])
    )
    df = similarity.executeStrategy()
    assert df.shape[0] == df.shape[1]
    for index, row in df.iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(report, str(index))
            gene2 = get_gene_list(report, str(col))
            expected = len(set(gene1) & set(gene2)) / min(len(gene1), len(gene2))
            assert np.isclose(value, expected, atol=1e-9), f"{index}, {col} invalid"


def test_cosine_protein():
    """タンパク質のcosine値を計算する"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import InteractionSimilarity
    from src.util.similarity_strategy import Cosine
    from tests.util.similarity import get_gene_list

    N_TEST = 20
    report = load_replicateIDR_report().head(N_TEST)

    similarity = InteractionSimilarity()
    similarity.setStrategy(
        Cosine(report=report, label_method=lambda df: df["Accession"])
    )
    df = similarity.executeStrategy()
    assert df.shape[0] == df.shape[1]
    for index, row in df.iterrows():
        for col, value in row.items():
            gene1 = get_gene_list(report, str(index))
            gene2 = get_gene_list(report, str(col))
            expected = len(set(gene1) & set(gene2)) / np.sqrt(len(gene1) * len(gene2))
            assert np.isclose(value, expected, atol=1e-9), f"{index}, {col} invalid"


def test_keywordcosine():
    """report.txtとsimilarity_matrixのラベル (タンパク質) が一致する"""
    """Multiple Sequence Analysis Distanceのテスト"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import InteractionSimilarity
    from src.util.similarity_strategy import KeywordCosine

    similarity = InteractionSimilarity()
    similarity.setStrategy(KeywordCosine())

    df = similarity.executeStrategy()
    print(df.head())

    expected = load_replicateIDR_report()["Target label"].unique().tolist()
    assert df.shape[0] == df.shape[1]
    assert len(df.columns) == len(expected)
    for protein in expected:
        assert protein in df.columns
        assert protein in df.index
    data = df.to_numpy()
    assert (data == data.T).reshape(-1).all(), "対称行列である"


def test_protein_transform():
    """Datasetへの変換ができる"""
    from src.util.bedfile import load_replicateIDR_report
    from src.util.similarity_protein import ProteinSimilarity
    from src.util.similarity_strategy import MSA

    N_TEST = 60
    report = load_replicateIDR_report().head(N_TEST)

    similarity = ProteinSimilarity()
    similarity.setStrategy(MSA(report=report))

    df = similarity.executeStrategy(transform=True)
    print(df.head())
    assert df.shape[0] == df.shape[1]
    assert len(df.columns) == len(report)
    data = df.to_numpy()
    assert (data == data.T).reshape(-1).all(), "対称行列である"
