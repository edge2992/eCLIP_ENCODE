def test_dataset_from_str():
    from src.util.bedfile import load_replicateIDR_report
    from src.eclip.dataset import Dataset

    report = load_replicateIDR_report()
    for index, row in report.head(3).iterrows():
        dataset = Dataset(row["Dataset"])
        assert dataset.dataset == row["Dataset"]
        assert dataset.protein == row["Target label"]
        assert dataset.biosample == row["Biosample name"]


def test_dataset_from_row():
    from src.util.bedfile import load_replicateIDR_report
    from src.eclip.dataset import Dataset

    report = load_replicateIDR_report()
    for index, row in report.head(3).iterrows():
        dataset = Dataset(row)
        assert dataset.dataset == row["Dataset"]
        assert dataset.protein == row["Target label"]
        assert dataset.biosample == row["Biosample name"]


def test_dataset_keywords():
    from src.util.bedfile import load_replicateIDR_report
    from src.eclip.dataset import Dataset

    report = load_replicateIDR_report().iloc[5]
    dataset = Dataset(report)
    assert len(dataset.keywords) == 9


def test_dataset_genes():
    from src.util.bedfile import load_replicateIDR_report
    from src.eclip.dataset import Dataset

    report = load_replicateIDR_report().iloc[5]
    dataset = Dataset(report)
    assert len(dataset.genes) == 32
