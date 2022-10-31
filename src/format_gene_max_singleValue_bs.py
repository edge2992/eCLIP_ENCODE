# 遺伝子ごとにsingleValueがmaxのbindingsiteを抽出して整形する
# prerequirements
# - downloaded_eCLIP_all
# - annotate_eCLIP_all
import os
import sys
from enum import Enum
from typing import Dict
import pandas as pd
from joblib import Parallel, delayed

sys.path.append("/mnt/H/MYWORK/eCLIP_ENCODE")
from src.annotate_eCLIP_all import get_output_file_path
from src.util.bedfile import (
    COLUMN_BED_NARROW_PEAK,
    PROJECT_PATH,
    load_report,
    read_intersected_bed,
)

PROJECT_PATH = "/mnt/H/MYWORK/eCLIP_ENCODE"


def transform_attribute_to_dict(attribute: str) -> Dict[str, str]:
    """gtfのattributeをdictに変換する"""
    splited_attributes = [x.strip() for x in attribute.replace('"', "").split(";")[:-1]]
    return dict(map(lambda x: x.split(" "), splited_attributes))


# use this for format_gene_binding_sites
FormatStrategy = Enum("FormatStrategy", ["MAX"])


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


def get_formatted_file_path(row: pd.Series, how):
    """整形済みのファイルのパスを取得する"""
    return os.path.join(
        PROJECT_PATH,
        "annotated_data",
        "gene",
        how.name.lower(),
        row["Assay term name"],
        row["Target label"],
        row["Biosample name"].split()[0],  # adrenal gland, K562, HepG2
        row["Accession"] + ".bed",
    )


def format_bed(row, how):
    """annotatedのbedfileを遺伝子ごとに整形する"""
    input_file = get_output_file_path(row)
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
