# eCLIPのピークファイルにgtfファイルを使って、
# そのピークがある遺伝子の情報を追加する
import itertools
import os
import subprocess
import sys
from typing import Dict, List

import pandas as pd
from dotenv import load_dotenv
from joblib import Parallel, delayed

load_dotenv()
PROJECT_PATH = os.environ["PROJECT_PATH"]

sys.path.append(PROJECT_PATH)
from src.util.get_bed_path import get_file_path

SORTED_DIR = os.path.join(PROJECT_PATH, "data/bed/sorted")
JACCARD_DIR = os.path.join(PROJECT_PATH, "data/bed/jaccard")


def sort_bed(input_path: str, save_path: str) -> subprocess.Popen:
    """bedfileをsortする"""
    cmd = "sort -k1,1 -k2,2n {} > {}".format(input_path, save_path)
    return subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)


def jaccard_bed(file_a: str, file_b: str, save_path: str) -> subprocess.Popen:
    cmd = "bedtools jaccard -a {} -b {} > {}".format(file_a, file_b, save_path)
    return subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE)


def sort_all(
    bedfiles: List[str], sorted_bedfiles: List[str], skip_if_exists: bool = True
):
    subprocs: Dict[str, subprocess.Popen] = {}
    for bedfile, sorted_bedfile in zip(bedfiles, sorted_bedfiles):
        if skip_if_exists and os.path.exists(sorted_bedfile):
            continue
        subprocs[os.path.basename(bedfile)] = sort_bed(bedfile, sorted_bedfile)

    for key, proc in subprocs.items():
        proc.wait()
        if proc.returncode != 0:
            print(proc.stderr.read())  # type: ignore
        else:
            print("{} Done".format(key))


def jaccard_bed_all(
    sorted_bedfiles: List[str], save_dir: str = JACCARD_DIR, skip_if_exists: bool = True
):
    subprocs: Dict[str, subprocess.Popen] = {}
    sorted_bedfiles.sort()
    for bed_a, bed_b in itertools.combinations(sorted_bedfiles, 2):
        accession_a = os.path.basename(bed_a).split(".")[0]
        accession_b = os.path.basename(bed_b).split(".")[0]
        save_path = os.path.join(save_dir, "{}_{}.txt".format(accession_a, accession_b))
        if skip_if_exists and os.path.exists(save_path):
            continue

        subprocs["{}_{}".format(accession_a, accession_b)] = jaccard_bed(
            bed_a, bed_b, save_path
        )

    for key, proc in subprocs.items():
        proc.wait()
        if proc.returncode != 0:
            print(proc.stderr.read())  # type: ignore
        else:
            print("{} Done".format(key))


def convert_df_from_jaccard_output(
    report_dir: str, metrics: str = "jaccard"
) -> pd.DataFrame:
    def load_jaccard_output(path: str) -> pd.Series:
        df = pd.read_table(path).iloc[0, :]
        assert isinstance(df, pd.Series)
        return df

    if metrics not in [
        "intersection",
        "union-intersection",
        "jaccard",
        "n-intersections",
    ]:
        raise ValueError(
            "metrics must be in [intersection, union-intersection, jaccard, n-intersections]"
        )
    files_file = [
        f for f in os.listdir(report_dir) if os.path.isfile(os.path.join(report_dir, f))
    ]

    series = Parallel(n_jobs=5, verbose=3)(
        delayed(load_jaccard_output)(os.path.join(report_dir, f)) for f in files_file
    )

    result = []

    for f, seri in zip(files_file, series):  # type: ignore
        basename = f.split(".")[0]
        dataset1 = basename.split("_")[0]
        dataset2 = basename.split("_")[1]
        result.append(
            {"dataset1": dataset1, "dataset2": dataset2, metrics: seri[metrics]}
        )
        result.append(
            {"dataset1": dataset2, "dataset2": dataset1, metrics: seri[metrics]}
        )
    data = pd.DataFrame(result)
    matrix = data.pivot(index="dataset1", columns="dataset2", values=metrics)
    for i in range(matrix.shape[0]):
        matrix.iloc[i, i] = 1
    return matrix


def main():
    from src.util.bedfile import load_replicateIDR_report

    report = load_replicateIDR_report().sort_values("Dataset")
    bedfiles = report.apply(lambda row: get_file_path(row), axis=1).tolist()
    sorted_bedfiles = [
        os.path.join(SORTED_DIR, os.path.basename(bedfile)) for bedfile in bedfiles
    ]

    if not os.path.exists(SORTED_DIR):
        os.makedirs(SORTED_DIR)

    sort_all(bedfiles, sorted_bedfiles)

    if not os.path.exists(JACCARD_DIR):
        os.makedirs(JACCARD_DIR)

    jaccard_bed_all(sorted_bedfiles, JACCARD_DIR)
    data = convert_df_from_jaccard_output(JACCARD_DIR)
    print(data.head())


if __name__ == "__main__":
    main()
