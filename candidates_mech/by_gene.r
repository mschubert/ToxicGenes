library(dplyr)
sys = import('sys')
util = import('../candidates/util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
    opt('o', 'outfile', 'pdf', 'by_gene/CDKN1A.rds')
)

cfg = yaml::read_yaml(args$config)
cohorts = c("BRCA", "COAD", "READ", "LUAD", "LUSC", "HNSC", "PRAD", "SKCM", "OV")

td = lapply(cohorts, util$load_tcga, top=args$gene) %>%
    bind_rows() %>%
    mutate(p53_mut = ifelse(is.na(p53_mut), "p53_wt", "p53_mut")) %>%
    group_by(cohort, gene) %>%
        mutate(expr = expr / max(expr, na.rm=TRUE)) %>%
    ungroup()

brca_sub = td %>%
    filter(cohort == "BRCA") %>%
    mutate(cohort = case_when(
        subtype %in% c("BRCA.LumA", "BRCA.LumB") ~ "BRCA.LumAB",
        subtype %in% c("BRCA.Basal") ~ "BRCA.Basal",
        TRUE ~ "other"
    )) %>%
    filter(cohort != "other")

td2 = bind_rows(td, brca_sub)

td3 = td2 %>%
    mutate(p53_mut = "all") %>%
    bind_rows(td2)

saveRDS(td3, file=args$outfile)
