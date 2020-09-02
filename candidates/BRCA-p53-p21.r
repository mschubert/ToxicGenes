library(dplyr)
library(ggplot2)
library(plyranges)
io = import('io')
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

idmap = import('process/idmap')
p53_targets = c("CDKN1A", "CCND2", "GADD45A", "FAS", "BAX", "CDKN2A", "TP53",
                "CCND3", "PERP", "CCNG1", "SESN1", "APAF1", "MDM2", "SERPINB5", "PIDD1", "ZMAT3")
rnaseq = tcga$rna_seq("BRCA", trans="vst")
rownames(rnaseq) = idmap$gene(rownames(rnaseq), to="hgnc_symbol")
rnaseq = rnaseq[p53_targets,]
p53_activity = tibble(sample=colnames(rnaseq), activity=scale(c(colMeans(rnaseq))))

prog_score = io$load("speed_linear.RData")
prog_score = tibble(sample=rownames(prog_score), prog_p53=prog_score[,"p53"])

#mut = tcga$mutations() %>%
#    filter(Study == "BRCA", Hugo_Symbol == "TP53") %>%
mut = io$load("/data/p282396/data/tcga/TCGAbiolinks-downloader/snv_mutect2/TCGA-BRCA.RData") %>%
    filter(Hugo_Symbol == "TP53") %>%
    transmute(Sample = substr(Tumor_Sample_Barcode, 1, 16),
              Variant = Variant_Classification,
              Var = Variant_Type,
              AAChange = Amino_acids) %>%
    mutate(p53_mut = case_when(
        Variant == "Missense_Mutation" ~ "p53_snp",
        Variant %in% c("Frame_Shift_Del", "Nonsense_Mutation", "Splice_Site") ~ "p53_del",
        TRUE ~ "p53_wt"))

brca_subtypes = readr::read_tsv("../../prmt5/BRCA.547.PAM50.SigClust.Subtypes.txt") %>%
    transmute(sample = substr(Sample, 1, 16),
              PAM50 = ifelse(PAM50 %in% "Basal", PAM50, "LumAB+Normal+HER2"))

td = util$load_tcga("BRCA", top="CDKN1A") %>%
    select(-p53_mut) %>% # what was that again?
    filter(cancer_copies >= 2-0.15) %>% # only euploid or amp
    as_tibble() %>%
    inner_join(brca_subtypes) %>%
    left_join(mut %>% dplyr::rename(sample=Sample, p53_var=Variant)) %>%
    left_join(p53_activity) %>%
    left_join(prog_score) %>%
    mutate(p53_mut = ifelse(is.na(p53_mut), "p53_wt", p53_mut),
           p53_var = ifelse(! p53_mut %in% c("p53_snp", "p53_del"), "other/none", p53_var),
           p53_var = relevel(factor(p53_var), "other/none"))
table(td$p53_mut)
table(td$p53_var)
#rlm3$do_fit()

tcga$intersect(td$sample, cpg, along=1)

#meth_cpg = tcga$meth_cpg("BRCA")

mean_eup = td %>% filter(cancer_copies < 2+0.15) %>%
    group_by(gene, p53_mut) %>%
    summarize(mean_eup = mean(expr))

pdf("BRCA-p53-p21.pdf", 10, 5)
ggplot(td, aes(x=cancer_copies, y=expr)) +
    geom_abline(data=mean_eup, aes(slope=mean_eup/2, intercept=0), linetype="dashed", color="red") +
    geom_point(aes(color=meth_eup_scaled, shape=p53_var), size=2, alpha=0.5) +
    geom_smooth(method="lm", se=FALSE, color="red") +
    scale_color_distiller(palette="RdBu") +
    facet_grid(PAM50 ~ p53_mut) +
    ggtitle("CDKN1A expression")

ggplot(td, aes(x=cancer_copies, y=meth)) +
    geom_point(aes(shape=p53_var, alpha=purity), size=2) +
    geom_smooth(method="lm", se=FALSE, color="blue") +
    facet_grid(PAM50 ~ p53_mut) +
    ggtitle("methylation")

ggplot(td, aes(x=cancer_copies, y=activity)) +
    geom_point(aes(shape=p53_var, alpha=purity), size=2) +
    geom_smooth(method="lm", se=FALSE, color="blue") +
    facet_grid(PAM50 ~ p53_mut) +
    ggtitle("p53 signature activity (geometric mean Veronica's genes)")

ggplot(td, aes(x=cancer_copies, y=prog_p53)) +
    geom_point(aes(shape=p53_var, alpha=purity), size=2) +
    geom_smooth(method="lm", se=FALSE, color="blue") +
    facet_grid(PAM50 ~ p53_mut) +
    ggtitle("p53 signature activity (progeny method)")

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
