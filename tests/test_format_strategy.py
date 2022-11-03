def test_format_max():
    """format_maxのテスト"""
    from src.util.bedfile import load_report, read_intersected_bed
    from src.util.get_bed_path import get_annotated_file_path
    from src.util.bed_format_strategy import format_max

    report = load_report()
    for _, row in report.head(3).iterrows():
        intersected = read_intersected_bed(get_annotated_file_path(row))
        result = format_max(intersected)
        counting_attribute = intersected["gff_attribute"].value_counts()
        testing_attribute = counting_attribute[counting_attribute > 1]
        for attribute, value in testing_attribute.head(10).items():
            if attribute == ".":
                continue
            testing_peaks = intersected[intersected["gff_attribute"] == attribute]
            gene_peak = result[result["gff_attribute"] == attribute]

            assert len(testing_peaks) == value
            assert len(gene_peak) == 1
            assert (
                testing_peaks["singleValue"].max() == gene_peak["singleValue"].values[0]
            )
