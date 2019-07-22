library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
sys = import('sys')

args = sys$cmd$parse(
    opt('s', 'sets', 'genes|gene set', 'genes'),
    opt('o', 'outfile', 'rds', 'merge_genes.rds'))

orf = readxl::read_xlsx(sprintf("../orf/pan/%s.xlsx", args$sets)) %>%
    mutate(fit = "lm")

ccle = list(
    rlm = readxl::read_xlsx(sprintf("../ccle/pan_rlm/%s.xlsx", args$sets)),
    rank = readxl::read_xlsx(sprintf("../ccle/pan_rank/%s.xlsx", args$sets))
) %>% bind_rows(.id="fit")

tcga_naive = list(
    rlm = readxl::read_xlsx(sprintf("../tcga/naive/pan_rlm/%s.xlsx", args$sets)),
    rank = readxl::read_xlsx(sprintf("../tcga/naive/pan_rank/%s.xlsx", args$sets))
) %>% bind_rows(.id="fit")
tcga_pur = list(
    rlm = readxl::read_xlsx(sprintf("../tcga/pur/pan_rlm/%s.xlsx", args$sets)),
    rank = readxl::read_xlsx(sprintf("../tcga/pur/pan_rank/%s.xlsx", args$sets))
) %>% bind_rows(.id="fit")
tcga_puradj = list(
    rlm = readxl::read_xlsx(sprintf("../tcga/puradj/pan_rlm/%s.xlsx", args$sets)),
    rank = readxl::read_xlsx(sprintf("../tcga/puradj/pan_rank/%s.xlsx", args$sets))
) %>% bind_rows(.id="fit")

tcga = list(
    naive = tcga_naive,
    pur = tcga_pur,
    puradj = tcga_puradj
) %>% bind_rows(.id="adj")

dset = list(orf=orf, ccle=ccle, tcga=tcga) %>%
    bind_rows(.id="dset") %>%
    select(name, dset, fit, adj, statistic) %>%
    mutate(adj = ifelse(is.na(adj), "none", adj))

saveRDS(dset, file=args$outfile)
