library(dplyr)
sys = import('sys')
util = import('../candidates/util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
    opt('o', 'outfile', 'pdf', 'by_gene/CDKN1A.rds')
)

cfg = yaml::read_yaml(args$config)
cohorts = c(setdiff(cfg$cor_tissues, c("pan", "NSCLC", "COADREAD")), "LUSC", "HNSC", "COAD")

td = lapply(cohorts, util$load_tcga, top=args$gene) %>%
    bind_rows() %>%
    mutate(p53_mut = ifelse(is.na(p53_mut), "p53_wt", "p53_mut")) %>%
    group_by(cohort, gene) %>%
        mutate(expr = expr / max(expr, na.rm=TRUE)) %>%
    ungroup()

saveRDS(td, file=args$outfile)
