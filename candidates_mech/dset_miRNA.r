library(dplyr)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')
util2 = import('./util')

#' Load miRNA expression data for a given gene
load_mirnas = function(cohort, gene) {
    mirnas = tcga$mirna_seq(cohort)
    mirnas = mirnas[,!duplicated(colnames(mirnas))] %>%
        DESeq2::DESeqDataSetFromMatrix(colData=data.frame(sample=colnames(.)), design=~ 1) %>%
        DESeq2::estimateSizeFactors() %>%
        DESeq2::counts(normalized=TRUE)

    binding = readRDS("../data/genesets/miRTarBase_2017.rds") %>%
        stack() %>%
        mutate(ind = sub("-[1-9]p?$", "", ind)) %>%
        filter(values == gene)

    #TODO: check if -[1-3] is the same miRNA or not (keeping it in for now)
    re = mirnas[sub("-[0-9]$" , "", rownames(mirnas)) %in% binding$ind,,drop=FALSE]
    rownames(re) = make.names(rownames(re))
    re
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/miRNA.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort)

    mirna = tryCatch(error = util2$muffle,
        lapply(cohorts, load_mirnas, gene=args$gene) %>%
            narray::stack(along=2) %>%
            tcga$map_id("specimen") %>% t())

    td = td %>% filter(sample %in% rownames(mirna))
    mirna = mirna[td$sample,]
    dset = cbind(td, mirna, constant=1)

    pdf(args$plotfile, 24, 12)
    print(plt$text(sprintf("miRNA expression (%i)", util2$nc(mirna)), size=20))
    for (v in colnames(mirna))
        print(util2$plot_l2d(dset, v, from=0))
    dev.off()
})
