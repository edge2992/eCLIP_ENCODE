import os
from time import sleep
from typing import List, Union

import pandas as pd
from dotenv import load_dotenv

from src.util.bedfile import load_replicateIDR_report
from stringdb import get_homology, get_interaction_partners, get_string_ids

load_dotenv()
PROJECT_PATH = os.getenv("PROJECT_PATH", "~/.cache/eclip")
SAVE_DIR = os.path.join(PROJECT_PATH, "data", "stringdb")


def identify_protein(
    query_protein: List[str], species: int = 9606, save_path: Union[str, None] = None
) -> pd.DataFrame:
    if save_path is None:
        save_path = os.path.join(SAVE_DIR, "stringdb_protein_id.tsv")

    if os.path.exists(save_path):
        data = pd.read_table(save_path)
    else:
        print("download STRINGdb protein id")
        data = get_string_ids(query_protein, species=species)
        print("save to", save_path)
        data.to_csv(save_path, sep="\t", index=False)

    return data


def download_homology_metrics(
    query_protein: List[str],
    species: int = 9606,
    save_path: Union[str, None] = None,
) -> pd.DataFrame:
    if save_path is None:
        save_path = os.path.join(SAVE_DIR, "stringdb_homology.tsv")

    if os.path.exists(save_path):
        data = pd.read_table(save_path)
    else:
        print("download STRINGdb homology metrics")
        data = get_homology(query_protein, species=species)
        print("save to", save_path)
        data.to_csv(save_path, sep="\t", index=False)

    return data


def download_interaction_partners(
    query_protein: List[str],
    species: int = 9606,
    required_score: int = 400,
    save_path: Union[str, None] = None,
) -> pd.DataFrame:
    if save_path is None:
        save_path = os.path.join(
            SAVE_DIR,
            "stringdb_interaction_partners_required_{}.tsv".format(required_score),
        )

    if os.path.exists(save_path):
        data = pd.read_table(save_path)
    else:
        print("download STRINGdb interaction partners")
        data = get_interaction_partners(
            query_protein, species=species, required_score=required_score
        )
        print("save to", save_path)
        data.to_csv(save_path, sep="\t", index=False)

    return data


if __name__ == "__main__":
    report = load_replicateIDR_report()
    query_protein = list(report["Target label"].unique())
    data = identify_protein(query_protein)
    homology = download_homology_metrics(data["preferredName"].to_list())
    print("homology: ", homology.shape)
    for required_score in [0, 400, 700, 900]:
        partners = download_interaction_partners(
            data["preferredName"].to_list(), required_score=required_score
        )
        print("partners: ", partners.shape)
        tmp = partners[
            partners["preferredName_B"].isin(partners["preferredName_A"].unique())
        ]
        print(required_score, ": ", partners.shape, " -> ", tmp.shape)
        sleep(1)
