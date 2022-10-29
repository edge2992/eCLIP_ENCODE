import os
import time
import requests

ENCODE_SEARCH_API = "https://www.encodeproject.org/search/?"
ENCODE_REPORT_TSV_API = "https://www.encodeproject.org/report.tsv?"
ENCODE_URL = "https://www.encodeproject.org"


def download_file(url, dir="data"):
    filename = os.path.join(dir, url.split("/")[-1])
    if os.path.exists(filename):
        print("File exists : {}".format(filename))
    else:
        print("Downloading: {}".format(url))
        cmd_wget = "wget -bqcN -P {} {}".format(dir, url)
        os.system(cmd_wget)
        time.sleep(0.25)


def download_table(url, dir="data"):
    filename = os.path.join(dir, "report.tsv")
    print("Downloading: {} -> {}".format(url, filename))
    cmd_wget = "wget -bqN -O {} '{}'".format(filename, url)
    os.system(cmd_wget)
    time.sleep(0.25)


def main(params):
    headers = {"accept": "application/json"}
    url = (
        ENCODE_SEARCH_API
        + "&".join([f"{k}={v}" for k, v in params.items()])
        + "&format=json"
    )
    response = requests.get(url, headers=headers)
    search_results = response.json()

    print(search_results["notification"], search_results["total"])
    assert search_results["notification"] == "Success"

    for entry in search_results["@graph"]:
        save_dir = os.path.join(
            "data",
            entry["assay_term_name"],
            entry["target"]["label"],
            entry["biosample_ontology"]["term_name"],
        )
        assert entry["assembly"] == "GRCh38"
        download_file(ENCODE_URL + entry["href"], save_dir)

    table_url = ENCODE_REPORT_TSV_API + "&".join(
        [f"{k}={v}" for k, v in params.items()]
    )
    download_table(table_url)


if __name__ == "__main__":
    params = {
        "type": "File",
        "assay_title": "eCLIP",
        "assembly": "GRCh38",
        "status": "released",
        "file_format": "bigBed",
        "file_type": "bigBed+narrowPeak",
        # "target.label": "AQR",
        "limit": "all",
    }
    main(params)
