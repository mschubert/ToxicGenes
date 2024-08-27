library(modules)
library(igraph)
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

wgd = readr::read_tsv("TCGA_mastercalls.abs_tables_JSedit.fixed.txt") %>%
    transmute(sample = paste0(array, "A"),
              wgd = as.integer(`Genome doublings`))

#mut = tcga$mutations()
#net = OmnipathR::interaction_graph(OmnipathR::import_all_interactions())
#as_ids(neighbors(net, "CDKN1A", "total"))

reads = lapply(cohorts, tcga$rna_seq) %>%
    narray::stack(along=2, allow_overwrite=TRUE) %>%
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
    inner_join(wgd) %>%
    na.omit() %>%
    dplyr::rename(sf = sizeFactor) %>%
    mutate(covar = tcga$barcode2study(sample),
           cancer_copies = pmax(0, (copies-2) / purity + 2))

saveRDS(df, file=args$outfile)
