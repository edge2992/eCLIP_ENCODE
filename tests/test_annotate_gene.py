def test_intersect_file_length():
    """bedfileが全て処理されていることを
    入力と出力のファイルの行数を比較して確認する"""
    from src.util.bedfile import load_report, count_file_length
    from src.util.get_bed_path import get_file_path
    from src.annotate_eCLIP_all import get_annotated_file_path
    from joblib import Parallel, delayed
    import pandas as pd

    def compare_file_length(row: pd.Series):
        input_file = get_file_path(row)
        output_file = get_annotated_file_path(row)
        return count_file_length(input_file) <= count_file_length(output_file)

    # report = load_report()
    report = load_report().head(10)

    result = Parallel(n_jobs=5, verbose=3)(
        delayed(compare_file_length)(row) for _, row in report.iterrows()
    )
    assert all(result)  # type: ignore
