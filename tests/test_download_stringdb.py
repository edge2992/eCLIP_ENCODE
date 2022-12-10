def test_idmapping():
    from src.util.download.stringdb import ENCODEprotein2preferredName

    mapping = ENCODEprotein2preferredName()
    samelist = [key for key in mapping.keys() if key in mapping.values()]
    print(len(mapping), " -> ", len(samelist))
    for key in samelist:
        assert mapping[key] == key
