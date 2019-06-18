library(dplyr)
library(cowplot)
plt = import('plot')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'overview.rds'),
    opt('o', 'outfile', 'xls', 'fits_naive.xlsx'),
    opt('p', 'plotfile', 'pdf', 'fits_naive.pdf'))

expr = readRDS(args$infile)

pan = expr %>%
    group_by(`GENE SYMBOL`, `Construct IDs`) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, function(d) broom::tidy(lm(`LFC DMSO/ETP` ~ 1, data=d)))) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    filter(term == "(Intercept)") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

pancov = expr %>%
    group_by(`GENE SYMBOL`, `Construct IDs`) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, function(d) broom::tidy(lm(`LFC DMSO/ETP` ~ tissue, data=d)))) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    filter(term == "(Intercept)") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

tissue = expr %>%
    mutate(`LFC DMSO/ETP` = `LFC DMSO/ETP` + runif(nrow(expr)) * 0.01) %>%
    group_by(`GENE SYMBOL`, `Construct IDs`, tissue) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, function(d) broom::tidy(lm(`LFC DMSO/ETP` ~ 1, data=d)))) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    filter(term == "(Intercept)") %>%
    select(-term) %>%
    group_by(tissue) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    arrange(adj.p, p.value)

result = list(pan=pan, pancov=pancov, tissue=tissue)

pdf(args$plotfile)
for (r in result)
    plt$volcano(r)
dev.off()

writexl::write_xlsx(args$outfile)
