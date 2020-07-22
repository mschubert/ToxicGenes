library(dplyr)
library(ggplot2)
library(plyranges)
seq = import('seq')
tcga = import('data/tcga')
util = import('./util')
#rlm3 = import('../tcga/fit_rlm3')

muffle = function(e) { warning(e, immediate.=TRUE); NULL }
load_cpg = function(cohort, gene) {
    cpgs = tcga$meth_cpg(cohort, annot=TRUE)
    idx = cpgs@rowRanges %>%
        filter(Gene_Symbol == gene)
    mat = SummarizedExperiment::assay(cpgs)[names(idx),]
}
cpg = tryCatch(error = muffle, { # in case no meth data
    re = lapply("BRCA", load_cpg, gene="CDKN1A") %>%
        narray::stack(along=2) %>%
        tcga$filter(primary=TRUE, cancer=TRUE) %>%
        tcga$map_id("specimen") %>% t()
    re = narray::map(re, along=1, subsets=tcga$barcode2study(rownames(re)),
    function(x) if (!all(is.na(x)) && sd(x, na.rm=TRUE)>0.05) x
        else rep(NA, length(x)))
    re = re[,narray::map(re, along=1, function(x) !all(is.na(x)))]
})
promoter_wanding = tcga$cpg_gene() %>%
    select(-everything()) %>%
    mutate(cg = names(.))
transcript_annots = seq$gene_table() %>%
    filter(external_gene_name == "CDKN1A") %>%
    mutate(chromosome_name = paste0("chr", chromosome_name),
           strand = setNames(c("+","-"),c(1,-1))[as.character(strand)]) %>%
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
    cg = intersect(colnames(cpg), gr_ids$cg)
    if (length(cg) > 0)
        narray::map(cpg[,cg,drop=FALSE], along=2, function(x) mean(x, na.rm=TRUE))
}
cgs2 = lapply(cgs, summarize_cgs) %>% do.call(cbind, .)
if (!is.null(cgs2) && ncol(cgs2) > 0)
    cpg = cbind(cgs2, cpg)

mut = tcga$mutations() %>%
    filter(Study == "BRCA", Hugo_Symbol == "TP53") %>%
    transmute(Sample = Tumor_Sample_Barcode,
              Variant = Variant_Classification,
              Var = Variant_Type,
              AAChange = AAChange)

snp = mut %>% filter(Var == "SNP") %>% pull(Sample) %>% unique()
del = mut %>% filter(Var == "DEL") %>% pull(Sample) %>% unique()

brca_subtypes = readr::read_tsv("../../prmt5/BRCA.547.PAM50.SigClust.Subtypes.txt") %>%
    transmute(sample = substr(Sample, 1, 16),
              PAM50 = ifelse(PAM50 %in% c("Basal", "Her2"), PAM50, "LumAB+Normal"))

td = util$load_tcga("BRCA", top="CDKN1A") %>%
    select(-p53_mut) %>%
    mutate(p53_mut = case_when(
        sample %in% snp ~ "p53_snp",
        sample %in% del ~ "p53_del",
        TRUE ~ "p53_wt")) %>%
    filter(cancer_copies >= 2-0.15) %>%
    as_tibble() %>%
    inner_join(brca_subtypes)
table(td$p53_mut)
#rlm3$do_fit()

tcga$intersect(td$sample, cpg, along=1)

#meth_cpg = tcga$meth_cpg("BRCA")

mean_eup = td %>% filter(cancer_copies < 2+0.15) %>%
    group_by(gene, p53_mut) %>%
    summarize(mean_eup = mean(expr))

pdf("BRCA-p53-p21.pdf", 10, 8)
ggplot(td, aes(x=cancer_copies, y=expr)) +
    geom_abline(data=mean_eup, aes(slope=mean_eup/2, intercept=0), linetype="dashed", color="red") +
    geom_point(aes(color=meth_eup_scaled), size=2, alpha=0.5) +
    geom_smooth(method="lm", se=FALSE, color="red") +
    scale_color_distiller(palette="RdBu") +
    facet_grid(PAM50 ~ p53_mut)

ggplot(td, aes(x=cancer_copies, y=meth)) +
    geom_point(aes(alpha=purity), size=2) +
    geom_smooth(method="lm", se=FALSE, color="blue") +
    facet_grid(PAM50 ~ p53_mut) +
    ggtitle("methylation")

td2 = cbind(td, cpg)
for (cp in colnames(cpg)) {
    p = ggplot(td2, aes_string(x="cancer_copies", y=cp)) +
        geom_point(aes(alpha=purity), size=2) +
        geom_smooth(method="lm", se=FALSE, color="blue") +
        facet_grid(PAM50 ~ p53_mut) +
        ggtitle(cp)
    print(p)
}
dev.off()
