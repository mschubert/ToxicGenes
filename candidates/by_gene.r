library(dplyr)
library(ggplot2)
library(patchwork)
library(plyranges)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
tcga = import('data/tcga')
util = import('./util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
    opt('p', 'plotfile', 'pdf', 'by_gene/CDKN1A.pdf'))

cfg = yaml::read_yaml(args$config)
cohorts = setdiff(cfg$cor_tissues, "pan")

# p53: some drop issues, should fix
td = lapply(cohorts, util$load_tcga, top=c("TP53", args$gene)) %>%
    bind_rows() %>%
    filter(cohort %in% cohorts, gene == args$gene)
#abl = util$summary_tcga(assocs, td) # add assocs?

### exon expression ###
load_exon = function(cohort, gene) {
    exons = tcga$rna_exon(cohort, annot=TRUE)
    emat = DESeq2::DESeqDataSet(exons, design=~1) %>%
        estimateSizeFactors() %>%
        counts(normalized=TRUE)
    idx = exons@rowRanges %>%
        filter(external_gene_name == gene)
    mat = emat[names(idx),]
}
exons = lapply(cohorts, load_exon, gene=args$gene) %>%
    narray::stack(along=2) %>%
    tcga$map_id("specimen") %>% t()

### cpg methylation ###
load_cpg = function(cohort, gene) {
    cpgs = tcga$meth_cpg(cohort, annot=TRUE)
    idx = cpgs@rowRanges %>%
        filter(Gene_Symbol == gene)
    mat = SummarizedExperiment::assay(cpgs)[names(idx),]
}
cpg = lapply(cohorts, load_cpg, gene=args$gene) %>%
    narray::stack(along=2) %>%
    tcga$map_id("specimen") %>% t()

### assemble dataset ###
colnames(exons) = make.names(colnames(exons))
dset = narray::stack(list(exons, cpg), along=2)
tcga$intersect(td$sample, dset, along=1)
cols = c("purity", colnames(exons), "meth", colnames(cpg))
dset = cbind(td, dset)

### plot ###
plot_l2d = function(dset, variable, et=0.15) {
    ggplot(dset, aes(x=cancer_copies, y=expr)) +
        util$stat_loess2d(aes_string(fill=variable), se_size=TRUE) +
        geom_density2d(bins=20, color="#00000050") +
        geom_vline(xintercept=c(2-et,2+et), color="black", linetype="dotted") +
        geom_point(data = td %>% filter(!is.na(mut)),
                   aes(shape=mut, size=is.na(mut)), color="black", alpha=1) +
        facet_wrap(~ cohort, scales="free", nrow=1) +
        scale_fill_viridis_c(option="magma", direction=-1) +
        scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                           values=c(0, seq_along(levels(td$mut))[-1]),
                           labels=levels(td$mut)) +
        guides(fill=guide_legend(title="")) +
        labs(title = variable,
             y = "normalized read count")
}

pdf(args$plotfile, 16, 4)
for (v in cols)
    print(plot_l2d(dset, v))
dev.off()
