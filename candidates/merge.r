library(dplyr)
library(patchwork)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('s', 'sets', 'genes|gene set', 'genes'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('o', 'outfile', 'rds', 'merge/pan/genes.rds'),
    opt('p', 'plotfile', 'pdf', 'merge/pan/genes.pdf'),
    opt('x', 'xlsx', 'xlsx', 'merge/pan/genes.xlsx'))

orf = readxl::read_xlsx(sprintf("../orf/%s/%s.xlsx", args$tissue, args$sets)) %>%
    mutate(adj = "none", fit = "lm", cna = "oe")

ccle = tidyr::crossing(adj = "none",
                       fit = c("rlm", "rlm2", "rank"),
                       cna = c("amp", "del", "all")) %>%
    mutate(data = purrr::pmap(list(fit, cna), function(fit, cna) {
        fname = sprintf("../ccle/%s_%s/%s.xlsx", args$tissue, fit, args$sets)
        message(fname)
        readxl::read_xlsx(fname, sheet = cna)
    })) %>%
    tidyr::unnest()

tcga = tidyr::crossing(adj = c("naive", "pur", "puradj"),
                       fit = c("rlm", "rlm2", "rank"),
                       cna = c("amp", "del", "all")) %>%
    mutate(data = purrr::pmap(list(adj, fit, cna), function(adj, fit, cna) {
        fname = sprintf("../tcga/%s/%s_%s/%s.xlsx", adj, args$tissue, fit, args$sets)
        message(fname)
        readxl::read_xlsx(fname, sheet = cna)
    })) %>%
    tidyr::unnest()

dset = list(orf=orf, ccle=ccle, tcga=tcga) %>%
    bind_rows(.id="dset") %>%
    select(name, dset, cna, fit, adj, estimate, statistic, adj.p) %>%
    group_by(dset, cna, fit, adj) %>%
        mutate(pctile = 100 * (1-rank(statistic)/n())) %>%
    ungroup() %>%
    mutate(dset = relevel(factor(dset), "orf"))

saveRDS(dset, file=args$outfile)

wide = dset %>%
    filter(adj %in% c("none", "puradj")) %>%
    transmute(name=name, subs=paste("fdr", dset, cna, sep="_"), sfdr=sign(statistic) * adj.p) %>%
    distinct(name, subs, .keep_all=TRUE) %>% # 19 dups pan gene
    tidyr::spread(subs, sfdr) %>%
    select(name, fdr_orf_oe, everything())
score = order(rowSums(1 - abs(wide[,-1]), na.rm=TRUE))
wide = wide[rev(score),]
writexl::write_xlsx(wide, args$xlsx)

venn = . %>% filter(estimate < 0, adj.p < 0.1) %>% select(name, dset) %>% distinct() %>% unstack()
pdf(args$plotfile)
plt$venn(venn(dset %>% filter(cna %in% c("oe", "all")))) + ggtitle("ALL compensated, 10% FDR")
plt$venn(venn(dset %>% filter(cna %in% c("oe", "amp")))) + ggtitle("AMP compensated, 10% FDR")
plt$venn(venn(dset %>% filter(cna %in% c("oe", "del")))) + ggtitle("DEL compensated, 10% FDR")
dev.off()
