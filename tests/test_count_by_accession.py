def test_convert_nested_dict_to_df():
    from src.util.bedfile import load_report, count_file_length
    from src.util.get_bed_path import get_file_path
    from src.plot.util.count_by_accession import (
        convert_gene_biosample_accession_count_dict_to_df,
        create_gene_biosample_accession_count_dict,
    )

    TEST_DATASET_NUM = 10
    report = load_report()
    dataset_list = report["Dataset"].unique()[:TEST_DATASET_NUM]
    small_report = report[report["Dataset"].isin(dataset_list)]
    nested_dict = create_gene_biosample_accession_count_dict(
        small_report,
        lambda row: count_file_length(get_file_path(row)),
    )
    df = convert_gene_biosample_accession_count_dict_to_df(nested_dict, small_report)
    assert df.shape == (TEST_DATASET_NUM, 3)
