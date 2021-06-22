library(dplyr)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
tcga = import('data/tcga')
idmap = import('process/idmap')
util2 = import('./util')

tf_cohort = function(cohort, tfs) {
    expr = tcga$rna_seq(cohort, trans="vst")
    rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
    expr[intersect(tfs, rownames(expr)),, drop=FALSE]
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/gsva-dorothea+tfexp.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort) %>% setdiff(c("BRCA.LumAB", "BRCA.Basal"))

    tfs = dorothea::dorothea_hs %>%
        filter(target == args$gene) %>%
        arrange(confidence) %>%
        mutate(label = sprintf("%s.%s_%s", tf, tolower(confidence),
                               setNames(c("act", "inh"), c("1", "-1"))[as.character(mor)]))

    scores = lapply(cohorts, tf_cohort, tfs=tfs$tf) %>%
        narray::stack(along=2) %>% t()
    scores = scores[,tfs$tf]
    colnames(scores) = tfs$label

    td = td %>% filter(sample %in% rownames(scores)) %>%
        mutate(cohort = ifelse(cohort %in% c("COAD", "READ"), "COADREAD", cohort),
               cohort = ifelse(cohort %in% c("LUAD", "LUSC"), "NSCLC", cohort))
    scores = scores[td$sample,]
    dset = cbind(td, scores, constant=1)

    pdf(args$plotfile, 28, 12)
    print(plt$text(sprintf("TF expression (%i)", util2$nc(scores)), size=20))
    for (v in colnames(scores))
        print(util2$plot_l2d(dset, v))
    dev.off()
})
