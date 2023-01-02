def test_count_gene_nunique():
    from src.util.bedfile import count_gene_nunique, load_report, read_annotated_bed
    from src.util.get_bed_path import get_formatted_file_path
    from src.util.bed_format_strategy import FormatStrategy

    report = load_report().head()
    for _, row in report.iterrows():
        filepath = get_formatted_file_path(row, FormatStrategy.MAX)
        df = read_annotated_bed(filepath)
        df = df[df["gene_id"].notnull()]
        num = count_gene_nunique(filepath)
        assert len(df) == num


def test_load_replicateIDR_drop_duplicate():
    from src.util.bedfile import load_replicateIDR_report

    report = load_replicateIDR_report(drop_duplicates=True)
    values = report["Date created"].map(lambda x: x.split("-")[0]).value_counts()
    print(values)
    assert values["2018"] == 218
    assert values["2022"] == 15
    assert values["2019"] == 2
