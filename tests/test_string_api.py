def test_get_network():
    from stringdb import get_network

    protein_list = "p35"
    result = get_network(protein_list)
    print(result)
    print(result[["preferredName_A", "preferredName_B", "score"]])


def test_get_network_image():
    from stringdb import get_network_image

    protein_list = "p35"
    result = get_network_image(protein_list)
    filepath = "tests/data/sample.png"
    with open(filepath, "rb") as fs:
        expected = fs.read()
    assert result == expected


def test_get_string_ids():
    from stringdb import get_string_ids

    protein_identifiers = "p35"
    result = get_string_ids(protein_identifiers)
    assert result.shape[0] == 1


def test_get_string_ids_list():
    from stringdb import get_string_ids

    protein_identifiers = ["p35", "cdk2"]
    result = get_string_ids(protein_identifiers)
    print(result[["queryIndex", "taxonName", "preferredName", "stringId"]])
    assert result.shape[0] == 2


def test_get_interaction_partners():
    from stringdb import get_interaction_partners

    protein_identifiers = ["p35", "cdk2"]
    result = get_interaction_partners(protein_identifiers)
    assert result.stringId_A.nunique() == len(protein_identifiers)


def test_get_homology():
    from stringdb import get_homology

    protein_identifiers = ["cdk1", "cdk2"]
    result = get_homology(protein_identifiers)
    assert result.shape[0] == 4
