def test_compare():
    from src.eclip import Dataset, Compare

    dataset1 = Dataset("/experiments/ENCSR384MWO/")
    dataset2 = Dataset("/experiments/ENCSR406OOZ/")
    compare = Compare([dataset1, dataset2])
    assert len(dataset1.genes) == 1428
    assert len(dataset2.genes) == 1590
    assert len(compare.keyword_intersection) == 6
    assert len(compare.gene_intersection) == 215
    assert (
        str(compare)
        == "/experiments/ENCSR384MWO/ (CSTF2, HepG2) & /experiments/ENCSR406OOZ/ (SUB1, HepG2)"
    )
    # dataset1: /experiments/ENCSR384MWO/
    # dataset2: /experiments/ENCSR406OOZ/
    # protein1: CSTF2
    # protein2: SUB1
    # gene1: 1428
    # gene2: 1590
    # keyword1: 11
    # keyword2: 12
    # keyword intersection: ['3D-structure', 'Reference proteome', 'Ubl conjugation', 'Isopeptide bond', 'Nucleus', 'Phosphoprotein']
    # intersection: 215
