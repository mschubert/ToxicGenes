library(dplyr)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
tcga = import('data/tcga')
idmap = import('process/idmap')
util = import('../candidates/util')
util2 = import('./util')

gsva_cohort = function(cohort, sets) {
    expr = tcga$rna_seq(cohort, trans="vst")
    rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
    GSVA::gsva(expr, sets, parallel.sz=0)
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/gsva-hallmarks.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort) %>% setdiff(c("BRCA.LumAB", "BRCA.Basal"))

    sets = gset$get_human("MSigDB_Hallmark_2020")
    scores = lapply(cohorts, gsva_cohort, sets=sets) %>%
        narray::stack(along=2) %>% t()
    colnames(scores) = make.names(colnames(scores))

    td = td %>% filter(sample %in% rownames(scores))
    scores = scores[td$sample,]
    dset = cbind(td, scores, constant=1)

    pdf(args$plotfile, 28, 12)
    print(plt$text(sprintf("MSigDB Hallmarks (%i)", util2$nc(scores)), size=20))
    for (v in colnames(scores))
        print(util2$plot_l2d(dset, v))
    dev.off()
})
