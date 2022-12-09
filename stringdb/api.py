from typing import List, Union
import requests
import os
import io
import pandas as pd

# STRING_API_BASE = "https://string-db.org/api"
STRING_API_BASE = "https://version-11-5.string-db.org"
CALLER_IDENTITY = "https://github.com/edge2992"


def request_url(endpoint: str, output_format="tsv") -> str:
    url = os.path.join(STRING_API_BASE, output_format, endpoint)
    return url


def handle_results(results: requests.Response):
    if results.ok:
        data = results.content.decode("utf-8")
        return pd.read_table(io.StringIO(data))
    else:
        raise ValueError(results.reason)


def get_network(
    identifiers: str,
    species: int = 9606,
    required_score: int = 400,
    caller_identity: str = CALLER_IDENTITY,
):
    params = {
        "identifiers": identifiers,
        "echo-query": 1,
        "species": species,
        "required_score": required_score,
        "caller_identity": caller_identity,
    }
    url = request_url("network")
    response = requests.post(url, data=params)
    return handle_results(response)


def get_homology(
    identifiers: Union[List, str],
    species: int = 9606,
    caller_identity: str = CALLER_IDENTITY,
):
    if isinstance(identifiers, list):
        identifiers = "\r".join(identifiers)
    params = {
        "identifiers": identifiers,
        "species": species,
        "caller_identity": caller_identity,
    }
    url = request_url("homology")
    response = requests.post(url, data=params)
    return handle_results(response)


def get_network_image(
    identifiers: str,
    species: int = 9606,
    add_white_nodes: int = 15,
    caller_identity: str = CALLER_IDENTITY,
):
    params = {
        "identifiers": identifiers,
        "echo-query": 1,
        "species": species,
        "add_white_nodes": add_white_nodes,
        "caller_identity": caller_identity,
    }
    url = request_url("network", output_format="image")
    response = requests.post(url, data=params)
    return response.content


def get_string_ids(
    identifiers: Union[List, str],
    species=9606,
    echo_query=1,
    caller_identity=CALLER_IDENTITY,
):
    if isinstance(identifiers, list):
        identifiers = "\r".join(identifiers)

    params = {
        "identifiers": identifiers,
        "species": species,
        "echo-query": echo_query,
        "caller_identity": caller_identity,
    }
    url = request_url("get_string_ids")
    response = requests.post(url, data=params)
    return handle_results(response)


def get_interaction_partners(
    identifiers: Union[str, List],
    species: int = 9606,
    required_score: int = 400,
    caller_identity: str = CALLER_IDENTITY,
):

    if isinstance(identifiers, list):
        identifiers = "\r".join(identifiers)

    params = {
        "identifiers": identifiers,
        "species": species,
        "required_score": required_score,
        "caller_identity": caller_identity,
    }
    url = request_url("interaction_partners")
    response = requests.post(url, data=params)
    return handle_results(response)
