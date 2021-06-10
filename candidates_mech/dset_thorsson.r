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
    cohorts = unique(td$cohort) %>% setdiff(c("BRCA.LumAB", "BRCA.Basal"))

    immune = load_thorsson()
    isubs = narray::mask(na.omit(setNames(immune$Immune.Subtype, rownames(immune)))) + 0
    measures = data.matrix(immune[!colnames(immune) %in% c("Immune.Subtype", "TCGA.Subtype")])
    td = td %>%
        mutate(cohort = ifelse(cohort %in% c("COAD", "READ"), "COADREAD", cohort),
               cohort = ifelse(cohort %in% c("LUAD", "LUSC"), "NSCLC", cohort)) %>%
        inner_join(data.frame(sample=rownames(immune), TCGA.Subtype=immune$TCGA.Subtype)) %>%
        filter(sample %in% rownames(isubs), sample %in% rownames(measures))

    isubs = isubs[td$sample,]
    measures = measures[td$sample,]
    dset = cbind(td, isubs, measures, constant=1)
    tsubs = split(dset, dset$cohort)

    pdf(args$plotfile, 28, 12)
    print(plt$text(sprintf("TCGA subtypes (%i)", length(tsubs), size=20)))
    for (cur in tsubs) {
        all = lapply(unique(cur$TCGA.Subtype),
                     function(x) cur %>% mutate(is_subtype=(TCGA.Subtype == x) + 0)) %>%
            bind_rows(.id="subtype")
        print(util2$plot_l2d(all, "is_subtype", facet=c("constant", "TCGA.Subtype")))
    }

    print(plt$text(sprintf("Immune measures (%i)", util2$nc(measures)), size=20))
    for (v in colnames(measures))
        print(util2$plot_l2d(dset, v))
    dev.off()

    print(plt$text(sprintf("Immune subtypes (%i)", util2$nc(isubs)), size=20))
    for (v in colnames(isubs))
        print(util2$plot_l2d(dset, v))
})
