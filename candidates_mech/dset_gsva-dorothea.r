library(dplyr)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
tcga = import('data/tcga')
idmap = import('process/idmap')
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
        opt('p', 'plotfile', 'pdf', 'CDKN1A/gsva-dorothea+target.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort) %>% setdiff(c("BRCA.LumAB", "BRCA.Basal"))

    sets = gset$get_human("DoRothEA")

    tfs = dorothea::dorothea_hs %>%
        mutate(confidence = factor(confidence, ordered=TRUE))
    p21 = tfs %>% filter(target == args$gene) %>% select(tf, conf_p21=confidence)
    sets = tfs %>%
        inner_join(p21) %>%
        mutate(label = sprintf("%s.%s_%s", tf, tolower(confidence),
                               setNames(c("act", "inh"), c("1", "-1"))[as.character(mor)])) %>%
        filter(confidence <= conf_p21, mor == 1) %>%
        select(target, label) %>%
        unstack() %>%
        gset$filter(min=3)

    scores = lapply(cohorts, gsva_cohort, sets=sets) %>%
        narray::stack(along=2) %>% t()

    td = td %>% filter(sample %in% rownames(scores)) %>%
        mutate(cohort = ifelse(cohort %in% c("COAD", "READ"), "COADREAD", cohort),
               cohort = ifelse(cohort %in% c("LUAD", "LUSC"), "NSCLC", cohort))
    scores = scores[td$sample,]
    dset = cbind(td, scores, constant=1)

    pdf(args$plotfile, 28, 12)
    print(plt$text(sprintf("DoRothEA TF activity (%i)", util2$nc(scores)), size=20))
    for (v in colnames(scores))
        print(util2$plot_l2d(dset, v))
    dev.off()
})
