import os
from typing import Dict, List, Union

import pandas as pd
from dotenv import load_dotenv

from src.util.bedfile import load_replicateIDR_report
from stringdb.api import get_homology, get_interaction_partners, get_string_ids

load_dotenv()
PROJECT_PATH = os.getenv("PROJECT_PATH", "~/.cache/eclip")
SAVE_DIR = os.path.join(PROJECT_PATH, "data", "stringdb")


def query_replicateIDR_protein() -> List:
    report = load_replicateIDR_report()
    return list(report["Target label"].unique())


def query_preferredName_replicateIDR_protein() -> List:
    data = identify_protein()
    return data["preferredName"].to_list()  # type: ignore


def ENCODEprotein2preferredName() -> Dict:
    from_protein = query_replicateIDR_protein()
    preferred_protein = query_preferredName_replicateIDR_protein()
    return dict(zip(from_protein, preferred_protein))


def identify_protein(
    species: int = 9606, save_path: Union[str, None] = None
) -> pd.DataFrame:
    if save_path is None:
        save_path = os.path.join(SAVE_DIR, "stringdb_protein_id.tsv")

    if os.path.exists(save_path):
        data = pd.read_table(save_path)
    else:
        query_protein = query_replicateIDR_protein()
        print("download STRINGdb protein id")
        data = get_string_ids(query_protein, species=species)
        print("save to", save_path)
        data.to_csv(save_path, sep="\t", index=False)

    return data


def download_homology_metrics(
    species: int = 9606,
    save_path: Union[str, None] = None,
) -> pd.DataFrame:
    if save_path is None:
        save_path = os.path.join(SAVE_DIR, "stringdb_homology.tsv")

    if os.path.exists(save_path):
        data = pd.read_table(save_path)
    else:
        print("download STRINGdb homology metrics")
        query_protein = query_preferredName_replicateIDR_protein()
        data = get_homology(query_protein, species=species)
        print("save to", save_path)
        data.to_csv(save_path, sep="\t", index=False)

    return data


def download_interaction_partners(
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
        query_protein = query_preferredName_replicateIDR_protein()
        data = get_interaction_partners(
            query_protein, species=species, required_score=required_score
        )
        print("save to", save_path)
        data.to_csv(save_path, sep="\t", index=False)

    return data
