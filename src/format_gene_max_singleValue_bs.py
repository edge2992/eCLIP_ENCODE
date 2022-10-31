# 遺伝子ごとにsingleValueがmaxのbindingsiteを抽出して整形する
# prerequirements
# - downloaded_eCLIP_all
# - annotate_eCLIP_all
import os
import sys
from typing import Dict
import pandas as pd
from joblib import Parallel, delayed
from dotenv import load_dotenv

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

sys.path.append(PROJECT_PATH)
from src.util.bedfile import (
    COLUMN_BED_NARROW_PEAK,
    load_report,
    read_intersected_bed,
)
from src.util.get_bed_path import get_annotated_file_path, get_formatted_file_path
from src.util.bed_format_strategy import FormatStrategy


def transform_attribute_to_dict(attribute: str) -> Dict[str, str]:
    """gtfのattributeをdictに変換する"""
    splited_attributes = [x.strip() for x in attribute.replace('"', "").split(";")[:-1]]
    return dict(map(lambda x: x.split(" "), splited_attributes))


def _format_gene_binding_sites(
    intersected: pd.DataFrame, how=FormatStrategy.MAX
) -> pd.DataFrame:
    """eCLIPのannotated peakを遺伝子ごとに整形する"""
    COLUMNS_INTENDED_ATTRIBUTES = ["gene_id", "gene_name", "gene_type"]
    formatted_peak: pd.DataFrame

    if how is FormatStrategy.MAX:
        formatted_peak = _format_max(intersected)
    else:
        raise ValueError("undefined format strategy")

    attribute = pd.DataFrame(
        formatted_peak["gff_attribute"].apply(transform_attribute_to_dict).to_list()  # type: ignore
    )
    return pd.concat(
        [
            formatted_peak[COLUMN_BED_NARROW_PEAK],
            attribute[COLUMNS_INTENDED_ATTRIBUTES],
        ],
        axis=1,
    )


def _format_max(intersected: pd.DataFrame):
    """singleValueがmaxのbindingsiteを抽出して整形する"""
    return (
        intersected.sort_values("singleValue", ascending=False)
        .drop_duplicates(["gff_seqname", "gff_start", "gff_end"], keep="first")
        .reset_index()
    )


def format_bed(row, how):
    """annotatedのbedfileを遺伝子ごとに整形する"""
    input_file = get_annotated_file_path(row)
    output_file = get_formatted_file_path(row, how)
    if os.path.exists(output_file):
        print("File exists: {}".format(output_file))
        return None
    else:
        dirpath = os.path.dirname(output_file)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        print("formatting: {}".format(os.path.basename(output_file)))
        intersected = read_intersected_bed(input_file)
        result = _format_gene_binding_sites(intersected, how)
        result.to_csv(output_file, sep="\t", index=False)


def main():
    report = load_report().head()
    how = FormatStrategy.MAX
    Parallel(n_jobs=5, verbose=3)(
        delayed(format_bed)(row, how) for _, row in report.iterrows()
    )


if __name__ == "__main__":
    main()
