library(dplyr)
library(GenomicRanges)
sys = import('sys')

args = sys$cmd$parse(opt('o', 'outfile', 'rds', 'exons.rds'))

ensembl = biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
exons = biomaRt::getBM(attributes=c("external_gene_name", "ensembl_exon_id", "exon_chrom_start",
                                    "exon_chrom_end", "chromosome_name", "strand"), mart=ensembl) %>%
    mutate(strand = setNames(c("+", "-"), c(1,-1))[as.character(strand)]) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)

saveRDS(exons, file=args$outfile)
