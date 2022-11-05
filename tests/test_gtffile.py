def test_gene_ids_gtf():
    from src.util.gtffile import gene_ids_gtf
    from src.util.bedfile import load_replicateIDR_report
    from src.plot.util.process_report import gene_ids_eCLIP

    genes = gene_ids_gtf()
    genes_eCLIP = gene_ids_eCLIP(load_replicateIDR_report())
    print("\ngff genes {} -> eCLIP genes {}".format(len(genes), len(genes_eCLIP)))
