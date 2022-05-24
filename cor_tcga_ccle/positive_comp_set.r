library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'orf_corr', 'xlsx', '../orf/fits_corrected.xlsx'),
    opt('c', 'ccle', 'xlsx', '../ccle/pan/stan-nb.xlsx'),
    opt('t', 'tcga_puradj', 'xlsx', '../tcga/pan/stan-nb_puradj.xlsx'),
    opt('o', 'outfile', 'tsv', 'positive_comp_set.tsv')
)

orf = readxl::read_xlsx(args$orf_corr)
ccle = readxl::read_xlsx(args$ccle)
tcga = readxl::read_xlsx(args$tcga_puradj)

both = inner_join(ccle %>% select(gene, est_ccle=estimate),
                  tcga %>% select(gene, est_tcga=estimate)) %>%
    mutate(hit = est_ccle < -0.3 & est_tcga < -0.3 & est_tcga + est_ccle < -0.8) %>%
    left_join(orf %>% select(gene=`GENE SYMBOL`, est_orf=estimate, stat_orf=statistic))

both$hit[is.na(both$hit)] = FALSE
table(both$hit)
"CDKN1A" %in% both$gene[both$hit]

write.table(both, file=args$outfile, sep="\t", row.names=FALSE, quote=FALSE)
