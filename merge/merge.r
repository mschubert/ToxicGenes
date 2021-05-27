library(dplyr)
library(ggplot2)
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('o', 'outfile', 'rds', 'pan.rds'),
    opt('p', 'plotfile', 'pdf', 'pan.pdf')
)

cfg = yaml::read_yaml(args$config)
#cnas = c("amp", "del", "all")
cnas = "amp"

orf = readxl::read_xlsx(sprintf("../orf/%s/%s.xlsx", args$tissue, "genes")) %>%
    mutate(adj = "none", fit = "lm", cna = "oe")

ccle = tidyr::crossing(adj = "none", fit = cfg$fits, cna = cnas) %>%
    mutate(data = purrr::pmap(list(fit, cna), function(fit, cna) {
        fname = sprintf("../ccle/%s/%s.xlsx", args$tissue, fit)
        message(fname)
        readxl::read_xlsx(fname) #, sheet = cna)
    })) %>%
    tidyr::unnest("data") %>%
    dplyr::rename(name = gene)

tcga = tidyr::crossing(adj = cfg$tcga_adj, fit = cfg$fits, cna = cnas) %>%
    mutate(data = purrr::pmap(list(adj, fit, cna), function(adj, fit, cna) {
        fname = sprintf("../tcga/%s/%s_%s.xlsx", args$tissue, fit, adj)
        message(fname)
        readxl::read_xlsx(fname) #, sheet = cna)
    })) %>%
    tidyr::unnest("data") %>%
    dplyr::rename(name = gene)

dset = list(orf=orf, ccle=ccle, tcga=tcga) %>%
    bind_rows(.id="dset") %>%
    select(name, dset, cna, fit, adj, estimate, eup_reads, statistic, p.value, adj.p) %>%
    group_by(dset, cna, fit, adj) %>%
        mutate(pctile = 100 * (1-rank(statistic)/n())) %>%
    ungroup() %>%
    mutate(dset = relevel(factor(dset), "orf"))

do_phist = function(dset, mod) {
    p = dset %>%
        filter(fit == mod) %>%
        mutate(dset = paste(dset, adj)) %>%
        ggplot(aes(x=p.value)) +
            geom_histogram(bins=50) +
            facet_grid(cna ~ dset)
}
plots = sapply(unique(dset$fit), do_phist, dset=dset, simplify=FALSE)
pdf(args$plotfile)
for (i in seq_along(plots))
    print(plots[[i]] + ggtitle(names(plots)[i]))
dev.off()

saveRDS(dset, file=args$outfile)
