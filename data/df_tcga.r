library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'df_tcga.rds')
)

cohorts = tcga$cohorts()
# excl x,y chroms

purity = tcga$purity() %>%
    filter(!is.na(estimate)) %>%
    transmute(sample = Sample,
              purity = estimate)

reads = lapply(cohorts, tcga$rna_seq) %>%
    narray::stack(along=2) %>%
    tcga$filter(cancer=TRUE, primary=TRUE)
rownames(reads) = idmap$gene(rownames(reads), to="hgnc_symbol")
reads = reads[rowMeans(reads) >= 10 & !is.na(rownames(reads)),]
eset = DESeq2::DESeqDataSetFromMatrix(reads, tibble(sample=colnames(reads)), ~1) %>%
    DESeq2::estimateSizeFactors() #todo: 16 duplicate hgnc_symbol

copies = lapply(cohorts, tcga$cna_genes) %>%
    narray::stack(along=2)
rownames(copies) = idmap$gene(rownames(copies), to="hgnc_symbol")
narray::intersect(reads, copies, along=1)
narray::intersect(reads, copies, along=2)
names(dimnames(reads)) = names(dimnames(copies)) = c("gene", "sample")

df = inner_join(reshape2::melt(reads, value.name="expr"),
                reshape2::melt(copies, value.name="copies")) %>%
    inner_join(purity) %>%
    inner_join(as.data.frame(SummarizedExperiment::colData(eset))) %>%
    na.omit() %>%
    dplyr::rename(sf = sizeFactor) %>%
    mutate(covar = tcga$barcode2study(sample),
           cancer_copies = pmax(0, (copies-2) / purity + 2),
           eup_dev = (copies - 2) / 2,
           eup_equiv = eup_dev + 1,
           eup_dev_cancer = (cancer_copies - 2) / 2,
           eup_equiv_cancer = eup_dev_cancer + 1)

saveRDS(df, file=args$outfile)
