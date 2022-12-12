import os
from time import sleep
from typing import List

from dotenv import load_dotenv

from src.util.bedfile import load_replicateIDR_report
from stringdb import get_network_image, get_string_ids

load_dotenv()
PROJECT_PATH = os.getenv("PROJECT_PATH", "~/.cache/eclip")
SAVE_DIR = os.path.join(PROJECT_PATH, "data", "stringdb", "network_images")


def download_network_images(save_dir: str, query_protein: List[str]):
    if os.path.exists(save_dir) is False:
        os.makedirs(save_dir)

    data = get_string_ids(query_protein)
    assert data.shape[0] == len(query_protein)

    for query, preferedName in zip(query_protein, data["preferredName"]):
        if query != preferedName:
            print(f"{query} -> {preferedName}")
        save_path = os.path.join(save_dir, f"{query}.png")
        if os.path.exists(save_path):
            print(f"skip {save_path}")
            continue
        image = get_network_image(preferedName)
        print(f"save to {save_path}")
        with open(save_path, "wb") as f:
            f.write(image)
        sleep(1)


if __name__ == "__main__":
    report = load_replicateIDR_report()
    query_protein = list(report["Target label"].unique())
    download_network_images(SAVE_DIR, query_protein)
