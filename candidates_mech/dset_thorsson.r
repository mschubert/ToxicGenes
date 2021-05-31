library(dplyr)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')
util2 = import('./util')

#' Load Thorsson et al. immune landscape data
load_thorsson = function() {
    immune_df = tcga$immune() %>%
        filter(cohort %in% cohorts) %>%
        select(-cohort, -OS, -`OS Time`, -PFI, -`PFI Time`) %>%
        as.data.frame()
    immune = immune_df[-1]
    rownames(immune) = paste0(immune_df$barcode, "-01A")
    colnames(immune) = make.names(colnames(immune))
    immune
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/thorsson.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort)

    immune = load_thorsson()
    isubs = narray::mask(na.omit(setNames(immune$Immune.Subtype, rownames(immune)))) + 0
    measures = data.matrix(immune[!colnames(immune) %in% c("Immune.Subtype", "TCGA.Subtype")])

    tcga$intersect(td$sample, isubs, measures, along=1)
    dset = cbind(td, isubs, measures, constant=1)

    pdf(args$plotfile, 24, 8)
    print(plt$text(sprintf("Immune measures (%i)", util2$nc(measures)), size=20))
    for (v in colnames(measures))
        print(util2$plot_l2d(dset, v))
    dev.off()

    print(plt$text(sprintf("Immune subtypes (%i)", util2$nc(isubs)), size=20))
    for (v in colnames(isubs))
        print(util2$plot_l2d(dset, v))
})
