library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
sys = import('sys')

args = sys$cmd$parse(
    opt('s', 'sets', 'genes|gene set', 'genes'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('o', 'outfile', 'rds', 'merge_genes.rds'))

orf = readxl::read_xlsx(sprintf("../orf/%s/%s.xlsx", args$tissue, args$sets)) %>%
    mutate(fit = "lm")

ccle = list(
    rlm = readxl::read_xlsx(sprintf("../ccle/%s_rlm/%s.xlsx", args$tissue, args$sets)),
    rank = readxl::read_xlsx(sprintf("../ccle/%s_rank/%s.xlsx", args$tissue, args$sets))
) %>% bind_rows(.id="fit")

tcga_naive = list(
    rlm = readxl::read_xlsx(sprintf("../tcga/naive/%s_rlm/%s.xlsx", args$tissue, args$sets)),
    rank = readxl::read_xlsx(sprintf("../tcga/naive/%s_rank/%s.xlsx", args$tissue, args$sets))
) %>% bind_rows(.id="fit")
tcga_pur = list(
    rlm = readxl::read_xlsx(sprintf("../tcga/pur/%s_rlm/%s.xlsx", args$tissue, args$sets)),
    rank = readxl::read_xlsx(sprintf("../tcga/pur/%s_rank/%s.xlsx", args$tissue, args$sets))
) %>% bind_rows(.id="fit")
tcga_puradj = list(
    rlm = readxl::read_xlsx(sprintf("../tcga/puradj/%s_rlm/%s.xlsx", args$tissue, args$sets)),
    rank = readxl::read_xlsx(sprintf("../tcga/puradj/%s_rank/%s.xlsx", args$tissue, args$sets))
) %>% bind_rows(.id="fit")

tcga = list(
    naive = tcga_naive,
    pur = tcga_pur,
    puradj = tcga_puradj
) %>% bind_rows(.id="adj")

dset = list(orf=orf, ccle=ccle, tcga=tcga) %>%
    bind_rows(.id="dset") %>%
    select(name, dset, fit, adj, statistic, adj.p) %>%
    mutate(adj = ifelse(is.na(adj), "none", adj)) %>%
    group_by(dset, fit, adj) %>%
        mutate(pctile = 100 * (1-rank(statistic)/n())) %>%
    ungroup() %>%
    mutate(dset = relevel(factor(dset), "orf"))

saveRDS(dset, file=args$outfile)
