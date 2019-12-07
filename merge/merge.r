library(dplyr)
library(patchwork)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('s', 'sets', 'genes|gene set', 'genes'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('o', 'outfile', 'rds', 'merge/pan/genes.rds'),
    opt('p', 'plotfile', 'pdf', 'merge/pan/genes.pdf'))

orf = readxl::read_xlsx(sprintf("../orf/%s/%s.xlsx", args$tissue, args$sets)) %>%
    mutate(adj = "none", fit = "lm", cna = "oe")

ccle = tidyr::crossing(adj = "none",
                       fit = c("rlm", "rlm3", "rank"),
                       cna = c("amp", "del", "all")) %>%
    mutate(data = purrr::pmap(list(fit, cna), function(fit, cna) {
        fname = sprintf("../ccle/%s_%s/%s.xlsx", args$tissue, fit, args$sets)
        message(fname)
        readxl::read_xlsx(fname, sheet = cna)
    })) %>%
    tidyr::unnest()

tcga = tidyr::crossing(adj = c("naive", "pur", "puradj"),
                       fit = c("rlm", "rlm3", "rank"),
                       cna = c("amp", "del", "all")) %>%
    mutate(data = purrr::pmap(list(adj, fit, cna), function(adj, fit, cna) {
        fname = sprintf("../tcga/%s/%s_%s/%s.xlsx", adj, args$tissue, fit, args$sets)
        message(fname)
        readxl::read_xlsx(fname, sheet = cna)
    })) %>%
    tidyr::unnest()

dset = list(orf=orf, ccle=ccle, tcga=tcga) %>%
    bind_rows(.id="dset") %>%
    select(name, dset, cna, fit, adj, estimate, eup_reads, rsq, statistic, adj.p) %>%
    group_by(dset, cna, fit, adj) %>%
        mutate(pctile = 100 * (1-rank(statistic)/n())) %>%
    ungroup() %>%
    mutate(dset = relevel(factor(dset), "orf"))

saveRDS(dset, file=args$outfile)
