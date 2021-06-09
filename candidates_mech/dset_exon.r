library(dplyr)
library(DESeq2)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')
util2 = import('./util')

#' Load exon expression data for a given gene
load_exon = function(cohort, gene) {
    exons = tcga$rna_exon(cohort, annot=TRUE)
    emat = DESeq2::DESeqDataSet(exons, design=~1) %>%
        estimateSizeFactors() %>%
        counts(normalized=TRUE)
    idx = exons@rowRanges %>%
        filter(external_gene_name == gene)
    mat = emat[names(idx),]
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/exon.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort)

    ### exon expression ###
    exons = tryCatch(error = util2$muffle, { # in case no exon data
        re = lapply(cohorts, load_exon, gene=args$gene) %>%
            narray::stack(along=2) %>%
            tcga$map_id("specimen") %>% t()
        colnames(re) = make.names(colnames(re))
        re
    })

    td = td %>% filter(sample %in% rownames(exons))
    exons = exons[td$sample,]
    dset = cbind(td, exons, constant=1)

    pdf(args$plotfile, 24, 12)
    print(plt$text(sprintf("Exon expression (%i)", util2$nc(exons)), size=20))
    for (v in colnames(exons))
        print(util2$plot_l2d(dset, v, from=0))
    dev.off()
})
