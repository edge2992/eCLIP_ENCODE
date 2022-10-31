# eCLIPのピークファイルにgtfファイルを使って、
# そのピークがある遺伝子の情報を追加する
import os
import sys
import subprocess
from typing import Union
from dotenv import load_dotenv

load_dotenv()
GENE_GTF_PATH = os.environ["GENE_GTF_PATH"]
PROJECT_PATH = os.environ["PROJECT_PATH"]

sys.path.append(PROJECT_PATH)
from src.util.bedfile import load_report
from src.util.get_bed_path import get_file_path, get_annotated_file_path


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


def main():
    if not os.path.exists(GENE_GTF_PATH):
        print("gtf file not found : {}".format(GENE_GTF_PATH))
        exit(1)
    report = load_report()
    subprocs = {}
    # 並列実行
    for _, row in report.iterrows():
        INPUT_BED_PATH = get_file_path(row)
        OUTPUT_BED_PATH = get_annotated_file_path(row)
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
