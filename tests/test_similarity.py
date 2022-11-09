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
    from src.util.similarity_protein import Similarity
    from src.util.similarity_strategy import MSA

    similarity = Similarity()
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
    from src.util.similarity_protein import Similarity
    from src.util.similarity_strategy import TAPE

    similarity = Similarity()
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
    from src.plot.util.process_report import label_protein_biosample, gene_ids_eCLIP
    from src.util.bedfile import read_annotated_bed
    from src.util.get_bed_path import get_formatted_file_path
    from src.util.bed_format_strategy import FormatStrategy
    from src.util.similarity_protein import Similarity
    from src.util.similarity_strategy import Lift
    import pandas as pd

    def get_gene_list(report: pd.DataFrame, label_pb: str):
        # TODO: 同じTarget label & Biosample nameで複数実験を行っている場合があるので、対応させる
        seri = report[
            (report["Target label"] == label_pb.split(maxsplit=1)[0])
            & (report["Biosample name"] == label_pb.split(maxsplit=1)[1])
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
    all_gene_count = len(gene_ids_eCLIP(report))

    similarity = Similarity()
    similarity.setStrategy(Lift(report=report))
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
