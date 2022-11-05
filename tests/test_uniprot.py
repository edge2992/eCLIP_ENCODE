def test_idmapping():
    from src.util.uniprot import idmapping
    from src.util.bedfile import load_replicateIDR_report

    eclip = load_replicateIDR_report()["Target label"].unique().tolist()
    d = idmapping()
    for protein in eclip:
        assert protein in d.keys()
