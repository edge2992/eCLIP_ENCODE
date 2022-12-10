import os
from time import sleep

from dotenv import load_dotenv

from src.util.bedfile import load_replicateIDR_report
from src.util.download.stringdb import (
    download_homology_metrics,
    download_interaction_partners,
    identify_protein,
)

load_dotenv()
PROJECT_PATH = os.getenv("PROJECT_PATH", "~/.cache/eclip")
SAVE_DIR = os.path.join(PROJECT_PATH, "data", "stringdb")


if __name__ == "__main__":
    report = load_replicateIDR_report()
    query_protein = list(report["Target label"].unique())
    data = identify_protein()
    homology = download_homology_metrics()
    print("homology: ", homology.shape)
    for required_score in [0, 400, 700, 900]:
        partners = download_interaction_partners(required_score=required_score)
        print("partners: ", partners.shape)
        tmp = partners[
            partners["preferredName_B"].isin(partners["preferredName_A"].unique())
        ]
        print(required_score, ": ", partners.shape, " -> ", tmp.shape)
        sleep(1)
