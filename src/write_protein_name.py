# タンパク質の名前を書き込む
import os
from typing import List
from src.util.bedfile import load_report
from dotenv import load_dotenv

load_dotenv()

PROJECT_PATH = os.environ["PROJECT_PATH"]


def write_proteins_name(proteins: List[str], filepath: str):
    with open(filepath, "w") as f:
        for protein in proteins:
            f.write(protein + "\n")


def main():
    proteins = list(load_report()["Target label"].unique())
    PROTEIN_PATH = os.path.join(PROJECT_PATH, "data", "uniprot", "proteins.txt")
    write_proteins_name(proteins, PROTEIN_PATH)


if __name__ == "__main__":
    main()
