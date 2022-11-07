# clustal omegaで距離行列を計算する
import os
from dotenv import load_dotenv

load_dotenv()
CLUSTALO_DISTMAT_PATH = os.environ["CLUSTALO_DISTMAT_PATH"]
CLUSTALO_CLUSTER_NUM_PATH = os.environ["CLUSTALO_CLUSTER_NUM_PATH"]


def run_clustalo(input_file: str):
    if CLUSTALO_DISTMAT_PATH == "":
        print("CLUSTALO_DISTMAT_PATH is not set")
    elif CLUSTALO_CLUSTER_NUM_PATH == "":
        print("CLUSTALO_CLUSTER_NUM_PATH is not set")
    elif os.path.exists(CLUSTALO_DISTMAT_PATH):
        print("File exists: {}".format(CLUSTALO_DISTMAT_PATH))
    else:
        dirpath = os.path.dirname(CLUSTALO_DISTMAT_PATH)
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        print("Running clustalo: {}".format(input_file))
        cmd_clustalo = (
            "clustalo --MAC-RAM 32000 --threads 10 --full --verbose "
            + "--outfmt clustal --resno --output-order input-order "
            + "-t Protein -i {} -o {} --distmat-out {} --force".format(
                input_file, CLUSTALO_CLUSTER_NUM_PATH, CLUSTALO_DISTMAT_PATH
            )
        )
        os.system(cmd_clustalo)


if __name__ == "__main__":
    INPUT_FASTA = "/mnt/H/MYWORK/eCLIP_ENCODE/data/uniprot/reviewed.fasta"
    run_clustalo(INPUT_FASTA)
