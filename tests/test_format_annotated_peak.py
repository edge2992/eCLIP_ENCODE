def test_format_peak_by_singleValue_max():
    """遺伝子ごとにsingleValue maxの値をその遺伝子のスコアとする"""
    from src.format_gene_max_singleValue_bs import (
        _format_gene_binding_sites,
        FormatStrategy,
    )
    from src.util.get_bed_path import get_annotated_file_path
    from src.util.bedfile import load_report, read_intersected_bed

    T_NUM = 5
    report = load_report().head(T_NUM)
    for _, row in report.iterrows():
        input_file = get_annotated_file_path(row)
        intersected = read_intersected_bed(input_file)
        result = _format_gene_binding_sites(intersected, FormatStrategy.MAX)

        print("{} {} -> {}".format(row["Accession"], intersected.shape, result.shape))
        print("gene_id ", len(result["gene_id"].unique()))
        print("gene_name ", len(result["gene_name"].unique()))
        # 異なるgene_nameが同じgene_idを持つ場合がある
        assert len(result["gene_id"].unique()) == result.shape[0]
