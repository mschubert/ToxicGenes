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

venn = function(title, df, comp=-0.5, fdr=0.01, r2=0.05) {
    df2 = df %>%
        filter(estimate < comp, adj.p < fdr, dset=="orf" | rsq > r2) %>%
        select(name, dset) %>% distinct()
    df3 = unstack(df2)
    df3 = df3[sapply(df3, length) > 0]
    p1 = plt$venn(df3) +
        ggtitle(sprintf("%s >= %i%% comp, %i%% FDR, %i%% R^2", title,
                        round(-comp*100), round(fdr*100), round(r2*100)))
    df4 = df2 %>%
        group_by(name) %>%
        mutate(n = n()) %>%
        ungroup() %>%
        arrange(-n) %>%
        filter(n >= 2) %>%
        mutate(name = factor(name, levels=rev(sort(unique(name)))))
    p2 = ggplot(df4, aes(x=factor(dset, levels=c("orf", "ccle", "tcga")), y=name)) +
        geom_tile(color="white", fill="black", size=0.5) +
        coord_fixed() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size=6, angle=90, hjust=1, vjust=0.5),
              axis.text.y = element_text(size=6))
    p1 + p2 + plot_layout(widths=c(5,1))
}
pdf(args$plotfile)
venn("amp+del", dset %>% filter(cna %in% c("oe", "all")))
venn("amp+del", dset %>% filter(cna %in% c("oe", "all")), r2=0.2)
venn("amp", dset %>% filter(cna %in% c("oe", "amp")))
venn("amp", dset %>% filter(cna %in% c("oe", "amp")), r2=0.2)
venn("del", dset %>% filter(cna %in% c("del")))
dev.off()
