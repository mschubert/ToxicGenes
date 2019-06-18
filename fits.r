library(dplyr)
library(cowplot)
plt = import('plot')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'overview.rds'),
    opt('f', 'field', 'response variable', 'LFC DMSO/ETP'),
    opt('o', 'outfile', 'xls', 'fits_naive.xlsx'),
    opt('p', 'plotfile', 'pdf', 'fits_naive.pdf'))

expr = readRDS(args$infile) %>%
    mutate(`LFC DMSO/ETP` = `LFC DMSO/ETP` + runif(nrow(.)) * 0.01)

y = rlang::sym(args$field)

pan = expr %>%
    group_by(`GENE SYMBOL`, `Construct IDs`) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, function(d) broom::tidy(lm(!! y ~ 1, data=d)))) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    filter(term == "(Intercept)") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

pancov = expr %>%
    group_by(`GENE SYMBOL`, `Construct IDs`) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, function(d) broom::tidy(lm(!! y ~ tissue, data=d)))) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    filter(term == "(Intercept)") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

tissue = expr %>%
    group_by(`GENE SYMBOL`, `Construct IDs`, tissue) %>%
    tidyr::nest() %>%
    mutate(result = purrr::map(data, function(d) broom::tidy(lm(!! y ~ 1, data=d)))) %>%
    select(-data) %>%
    tidyr::unnest() %>%
    filter(term == "(Intercept)") %>%
    select(-term) %>%
    group_by(tissue) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    ungroup() %>%
    arrange(adj.p, p.value)

result = list(pan=pan, pancov=pancov, tissue=tissue)

plots = lapply(result, function(res) {
    if ("tissue" %in% names(res)) {
        res = mutate(res, label = sprintf("%s - %s", `GENE SYMBOL`, tissue))
    } else {
        res = mutate(res, label = `GENE SYMBOL`)
    }

    if (all(res$estimate[rank(res$p.value) < 10] > 0))
        res$label[res$estimate > 0] = NA

    res %>%
        mutate(size=5) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", thresh=0.1, dir=-1) %>%
        plt$volcano(label_top=30, repel=TRUE)
})

pdf(args$plotfile)
for (i in seq_along(plots))
    print(plots[[i]] + ggtitle(names(plots)[i]))
dev.off()

writexl::write_xlsx(result, args$outfile)
