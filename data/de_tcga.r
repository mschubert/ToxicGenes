library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')
deseq = import('process/deseq')

normal_cancer_de = function(cohort) {
    message(cohort)
    reads = tcga$rna_seq(cohort) %>% tcga$filter(primary=TRUE)
    reads = reads[,!duplicated(colnames(reads))]
    reads = reads[,substr(colnames(reads), 14, 16) %in% c("01A", "11A")]
    cols = tcga$barcode2index(colnames(reads)) %>%
        mutate(Sample.Definition = factor(Sample.Definition,
            levels=c("Solid Tissue Normal", "Primary Solid Tumor")))
    if (any(cols$Sample.Definition == "Solid Tissue Normal"))
        DESeq2::DESeqDataSetFromMatrix(reads, cols, ~Sample.Definition) %>%
            DESeq2::estimateSizeFactors() %>%
            DESeq2::DESeq() %>%
            deseq$extract_result() %>%
            arrange(padj, pvalue)
}

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'de_tcga.rds')
)

cohorts = tcga$cohorts()
res = sapply(cohorts, normal_cancer_de, simplify=FALSE)
res = res[!sapply(res, is.null)]

saveRDS(res, file=args$outfile)
