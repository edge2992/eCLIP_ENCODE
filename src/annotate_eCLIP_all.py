# eCLIPのピークファイルにgtfファイルを使って、
# そのピークがある遺伝子の情報を追加する
import pandas as pd
import os
import sys
import subprocess
from typing import Union

sys.path.append("/mnt/H/MYWORK/eCLIP_ENCODE")
from src.plot.util.bedfile import PROJECT_PATH, get_file_path, load_report


GENE_GTF_PATH = "/home/edge2992/Resource/gencode.v24.annotation.gene.gtf"


def intersect_bed(input_file: str, output_file: str) -> Union[subprocess.Popen, None]:
    """intersectBedを使って、bedファイルをgtfファイルとintersectする"""
    dirpath = os.path.dirname(output_file)
    if os.path.exists(output_file):
        print("File exists :  {}".format(output_file))
        return None
    else:
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        print("intersecting: {}".format(os.path.basename(input_file)))
        cmd_intersect = "intersectBed -a {} -b {} -wa -s -loj > {}".format(
            input_file, GENE_GTF_PATH, output_file
        )
        return subprocess.Popen(cmd_intersect, shell=True, stderr=subprocess.PIPE)


def get_output_file_path(row: pd.Series):
    """report.tsvの行からファイルのパスを取得する"""
    return os.path.join(
        PROJECT_PATH,
        "annotated_data",
        "intersect",
        row["Assay term name"],
        row["Target label"],
        row["Biosample name"].split()[0],  # adrenal gland, K562, HepG2
        row["Accession"] + ".bed",
    )


def main():
    if not os.path.exists(GENE_GTF_PATH):
        print("gtf file not found : {}".format(GENE_GTF_PATH))
        exit(1)
    report = load_report()
    subprocs = {}
    # 並列実行
    for _, row in report.iterrows():
        INPUT_BED_PATH = get_file_path(row)
        OUTPUT_BED_PATH = get_output_file_path(row)
        proc = intersect_bed(INPUT_BED_PATH, OUTPUT_BED_PATH)
        if proc is not None:
            subprocs[row["Accession"]] = proc

    # 回収
    for key, proc in subprocs.items():
        proc.wait()
        if proc.returncode != 0:
            print(proc.stderr.read())
        else:
            print("{} Done".format(key))


if __name__ == "__main__":
    main()
