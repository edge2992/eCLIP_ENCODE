def test_convert_nested_dict_to_df():
    from src.util.bedfile import load_report, count_file_length
    from src.util.get_bed_path import get_file_path
    from src.plot.util.process_by_accession import (
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


def test_convert_accession_count_length_dict():
    """アッセイごとにファイル数 (ピーク数) をカウントする"""
    from src.util.bedfile import load_report, count_file_length
    from src.util.get_bed_path import get_file_path
    from src.plot.util.process_by_accession import _create_accession_count_dict

    expected = {
        "ENCFF431NBT": 38156,
        "ENCFF645PLG": 26254,
        "ENCFF913FPP": 8373,
        "ENCFF419MOS": 219460,
        "ENCFF394EQW": 81236,
        "ENCFF402AIE": 92058,
        "ENCFF384NHH": 9279,
        "ENCFF371SNO": 57572,
        "ENCFF210TQC": 201726,
        "ENCFF067JAD": 255064,
    }
    TEST_DATASET_NUM = 10
    report = load_report()
    result = _create_accession_count_dict(
        report.head(TEST_DATASET_NUM),
        lambda row: count_file_length(get_file_path(row)),
    )
    assert len(result) == TEST_DATASET_NUM
    for key, value in result.items():
        assert value == expected[key]


def test_convert_accession_count_gene_dict():
    """アッセイごとに遺伝子数をカウントする"""
    from src.util.bedfile import load_report
    from src.util.get_bed_path import get_formatted_file_path
    from src.plot.util.process_by_accession import _create_accession_count_dict
    from src.util.bedfile import count_gene_nunique
    from src.util.bed_format_strategy import FormatStrategy

    expected = {
        "ENCFF431NBT": 7949,
        "ENCFF645PLG": 5560,
        "ENCFF913FPP": 3326,
        "ENCFF419MOS": 14125,
        "ENCFF394EQW": 9163,
        "ENCFF402AIE": 10505,
        "ENCFF384NHH": 3702,
        "ENCFF371SNO": 9012,
        "ENCFF210TQC": 14253,
        "ENCFF067JAD": 13562,
    }
    TEST_DATASET_NUM = 10
    report = load_report()
    result = _create_accession_count_dict(
        report.head(TEST_DATASET_NUM),
        lambda row: count_gene_nunique(
            get_formatted_file_path(row, FormatStrategy.MAX)
        ),
    )
    assert len(result) == TEST_DATASET_NUM
    for key, value in result.items():
        assert value == expected[key]


def test_convert_accession_gene_dict():
    """アッセイごとに遺伝子のリストを作成する"""
    from src.util.bedfile import load_report
    from src.plot.util.process_report import get_gene_ids

    expected = {
        "ENCFF431NBT": 7949,
        "ENCFF645PLG": 5560,
        "ENCFF913FPP": 3326,
        "ENCFF419MOS": 14125,
        "ENCFF394EQW": 9163,
        "ENCFF402AIE": 10505,
        "ENCFF384NHH": 3702,
        "ENCFF371SNO": 9012,
        "ENCFF210TQC": 14253,
        "ENCFF067JAD": 13562,
    }

    TEST_DATASET_NUM = 10
    report = load_report()
    result = get_gene_ids(
        report.head(TEST_DATASET_NUM),
    )
    assert len(result) == TEST_DATASET_NUM
    for key, value in result.items():
        assert len(value) == expected[key]
