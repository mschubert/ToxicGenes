library(dplyr)
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')
util2 = import('./util')

get_ensembl_rest_info = function(gene_id) {
    query = sprintf("https://rest.ensembl.org/lookup/id/%s?expand=1", gene_id)
    r = httr::GET(query, httr::content_type("application/json"))
    httr::stop_for_status(r)
    res = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))

    res$Transcript %>% select(ensembl_transcript_id=id, Exon) %>% tidyr::unnest(Exon) %>%
        lapply(unlist) %>% do.call(data.frame, .) %>% as_tibble()
}

plot_gene_annot = function(gene_name, cpg) {
    trs = seq$coords$transcript() %>%
        filter(external_gene_name == gene_name)
    strand = setNames(c("-", "+"), c("-1","1"))[as.character(trs$strand[1])]
    adds = get_ensembl_rest_info(trs$ensembl_gene_id[1]) %>%
        select(label=ensembl_transcript_id, exon_id=id, start, end) %>%
        mutate(type="transcript")
    trs2 = trs %>%
        transmute(label=ensembl_transcript_id, start=transcript_start, end=transcript_end)

    probes = grep("^cg", colnames(cpg), value=TRUE)
    cgs = tcga$meth_cg2gene() %>%
        filter(probeID %in% probes) %>%
        transmute(label=probeID, start=CpG_beg, end=CpG_end) %>%
        left_join(tibble(label=colnames(cpg), cg_sd=apply(cpg, 2, sd)))

    both = bind_rows(list(transcript=trs2, cgs=cgs), .id="type")

    ggplot(both, aes(x=start, y=label, color=type)) +
        geom_segment(aes(xend=end, yend=label), size=1) +
        geom_segment(data=adds, aes(xend=end, yend=label), size=4, alpha=0.5) +
        geom_point(aes(size=cg_sd)) +
        ggtitle(sprintf("%s (%s)", gene_name, strand))
}

meth_for_gene = function(cohorts, gene_name, gene_id="external_gene_name") {
    cg1 = tcga$meth_summary(cohorts, gene_id)[gene_name,,]
    map = tcga$meth_mapping(gene_id)$pgene
    probe_ids = map$probe_id[map[[gene_id]] == gene_name]
    cg2 = lapply(cohorts, function(c) {
        m = tcga$meth_cpg(c)
        m[intersect(rownames(m), probe_ids),,drop=FALSE]
    }) %>% narray::stack(along=2) %>% t()
    cg2 = cg2[,names(sort(-apply(cg2, 2, sd)))]
    re = narray::stack(cg1, cg2, along=2)
    re[,colSums(is.na(re)) < nrow(re)/2]
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/meth.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort) %>% setdiff(c("BRCA.LumAB", "BRCA.Basal"))
    cpg = tryCatch(meth_for_gene(cohorts, args$gene), error=util2$muffle)

    td = td %>% filter(sample %in% rownames(cpg)) %>%
        mutate(cohort = ifelse(cohort %in% c("COAD", "READ"), "COADREAD", cohort),
               cohort = ifelse(cohort %in% c("LUAD", "LUSC"), "NSCLC", cohort))
    cpg = cpg[td$sample,]
    dset = cbind(td, cpg, constant=1)

    pdf(args$plotfile, 28, 12)
    print(plt$text(sprintf("Methylation (%i CpG)", util2$nc(cpg)), size=20))
    print(plot_gene_annot(args$gene, cpg))
    if (sum(!is.na(dset$meth_eup_scaled)) > 0)
        print(util2$plot_l2d(dset, "meth_eup_scaled"))

    for (v in colnames(cpg))
        print(util2$plot_l2d(dset, v))
    dev.off()
})
