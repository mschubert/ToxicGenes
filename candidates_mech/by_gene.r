library(dplyr)
library(ggplot2)
library(patchwork)
library(plyranges)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')

#' Handle error as warning and return NULL
muffle = function(e) { warning(e, immediate.=TRUE); NULL }

#' Get number of cases in matrix
nc = function(mat) { re = ncol(mat); if (is.null(re)) 0 else re }

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
    opt('p', 'plotfile', 'pdf', 'by_gene/CDKN1A.pdf'))

cfg = yaml::read_yaml(args$config)
cohorts = c(setdiff(cfg$cor_tissues, "pan"), "LUSC", "HNSC", "COAD")

# p53: some drop issues, should fix
td = lapply(cohorts, util$load_tcga, top=c("TP53", args$gene)) %>%
    bind_rows() %>%
    filter(cohort %in% cohorts, gene == args$gene) %>%
    mutate(p53_mut = ifelse(is.na(p53_mut), "p53_wt", "p53_mut")) %>%
    group_by(cohort, gene) %>%
        mutate(expr = expr / max(expr, na.rm=TRUE)) %>%
    ungroup()

### rna isoforms ###
load_isoform = function(cohort, gene) {
}
#isos = lapply(cohorts, load_isoform, gene=args$gene)

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
exons = tryCatch(error = muffle, { # in case no exon data
    re = lapply(cohorts, load_exon, gene=args$gene) %>%
        narray::stack(along=2) %>%
        tcga$map_id("specimen") %>% t()
    colnames(re) = make.names(colnames(re))
    re
})

### miRNAs ###
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
mirna = tryCatch(error = muffle,
    lapply(cohorts, load_mirnas, gene=args$gene) %>%
        narray::stack(along=2) %>%
        tcga$map_id("specimen") %>% t())

### cpg methylation ###
load_cpg = function(cohort, gene) {
    cpgs = tcga$meth_cpg(cohort, annot=TRUE)
    idx = cpgs@rowRanges %>%
        filter(Gene_Symbol == gene)
    mat = SummarizedExperiment::assay(cpgs)[names(idx),]
}
cpg = tryCatch(error = muffle, # in case no meth data
    lapply(cohorts, load_cpg, gene=args$gene) %>%
        narray::stack(along=2) %>%
        tcga$map_id("specimen") %>% t())

promoter_wanding = tcga$cpg_gene() %>%
    transmute(cg = names(.))
transcript_annots = seq$gene_table() %>%
    filter(external_gene_name == args$gene) %>%
    mutate(chromosome_name = paste0("chr", chromosome_name),
           strand = setNames(c("+","-"),c(1,-1))[strand]) %>%
    GenomicRanges::makeGRangesFromDataFrame()
cgs = list(
    core = transcript_annots %>% anchor_5p() %>% mutate(width=400) %>% shift_upstream(200) %>%
        join_overlap_intersect(promoter_wanding),
    ext = transcript_annots %>% anchor_5p() %>% mutate(width=3000) %>% shift_upstream(1500) %>%
        join_overlap_intersect(promoter_wanding),
    body = transcript_annots %>% filter(width > 1500) %>% anchor_3p() %>% stretch(-1500) %>%
        join_overlap_intersect(promoter_wanding)
)
summarize_cgs = function(gr_ids) {
    cgs = intersect(colnames(cpg), gr_ids$cg)
    if (length(cgs) > 0)
        narray::map(cpg[,cgs], along=2, function(x) mean(x, na.rm=TRUE))
}
cgs2 = lapply(cgs, summarize_cgs) %>% do.call(cbind, .)
if (!is.null(cgs2) && ncol(cgs2) > 0)
    cpg = cbind(cgs2, cpg)

### cell types ###
immune_df = tcga$immune() %>%
    filter(cohort %in% cohorts) %>%
    select(-cohort, -`Immune Subtype`, -`TCGA Subtype`, -OS, -`OS Time`, -PFI, -`PFI Time`)
immune = data.matrix(immune_df[-1])
rownames(immune) = paste0(immune_df$barcode, "-01A")
colnames(immune) = make.names(colnames(immune))

### assemble dataset ###
dset = narray::stack(list(exons, cpg, mirna, immune), along=2)
tcga$intersect(td$sample, dset, along=1)
dset = cbind(td, dset)

### plot ###
plot_l2d = function(dset, variable, et=0.15, from=NA, to=NA) {
    if (all(na.omit(dset[[variable]]) >= 0))
        fill = scale_fill_viridis_c(option="magma", direction=-1, limits=c(from, to))
    else
        fill = scale_fill_gradientn(colours=rev(RColorBrewer::brewer.pal(11,"RdBu")),
                                limits=c(from, to))
    ggplot(dset, aes(x=cancer_copies, y=expr)) +
        util$stat_gam2d(aes_string(fill=variable, by="purity"), se_size=TRUE) +
        geom_density2d(bins=20, color="chartreuse4", size=0.7) +
        geom_vline(xintercept=c(2-et,2+et), color="springgreen4", linetype="dotted", size=1.5) +
        facet_grid(p53_mut ~ cohort, scales="free") +
        fill +
        scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                           values=c(0, seq_along(levels(td$mut))[-1]),
                           labels=levels(td$mut)) +
        guides(fill=guide_legend(title="")) +
        labs(title = variable,
             y = "normalized read count")
}

pdf(args$plotfile, 24, 8)
print(plot_l2d(dset, "purity"), from=0, to=1)
print(plot_l2d(dset, "expr", from=0))

print(plt$text(sprintf("Exon expression (%i)", nc(exons)), size=20))
for (v in colnames(exons))
    print(plot_l2d(dset, v, from=0))

print(plt$text(sprintf("miRNA expression (%i)", nc(mirna)), size=20))
for (v in colnames(mirna))
    print(plot_l2d(dset, v, from=0))

print(plt$text(sprintf("Methylation (%i CpG)", nc(cpg)), size=20))
if (sum(!is.na(dset$meth_eup_scaled)) > 0)
    print(plot_l2d(dset, "meth_eup_scaled"))

for (v in colnames(cpg))
    print(plot_l2d(dset, v, from=0, to=1))

print(plt$text(sprintf("Immune subtypes (%i)", nc(immune)), size=20))
for (v in colnames(immune))
    print(plot_l2d(dset, v))

dev.off()
