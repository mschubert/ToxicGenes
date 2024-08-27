library(modules)
library(dplyr)
library(ggplot2)
plt = import('plot')
sys = import('sys')

do_fit = function(expr, fml) {
    ffun = function(g) {
        expr %>%
            mutate(gene = `GENE SYMBOL` == g) %>%
            lm(fml, data=.) %>%
            broom::tidy() %>%
            mutate(size = sum(expr$`GENE SYMBOL` == g))
    }

    res = expr %>%
        select(`GENE SYMBOL`) %>%
        distinct() %>%
        mutate(result = clustermq::Q(ffun, `GENE SYMBOL`, n_jobs=10,
            export=list(expr=expr, fml=fml), pkgs="dplyr")) %>%
        tidyr::unnest() %>%
        filter(term == "geneTRUE") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'overview.rds'),
        opt('o', 'outfile', 'xls', 'fitsWGD.rds')
    )

    expr = readRDS(args$infile) %>%
        mutate(`LFC DMSO/ETP` = `LFC DMSO/ETP` + runif(nrow(.)) * 0.01)

    `WGD+` = c("BT474", "D458", "Kuramochi", "NB69", "SK-BR-3", "OVSAHO", "T47D", "WM266-4")
    `WGD-` = c("Meljuso", "TC32")
    # NA: D283, LnCaP, OVCAR4
    # not in data: BE2C, H2077, LAN-1, SK-NEP-1

    result = list(
        `panWGD+` = do_fit(expr[expr$cells %in% `WGD+`,], `LFC DMSO/ETP` ~ gene),
        `panWGD-` = do_fit(expr[expr$cells %in% `WGD-`,], `LFC DMSO/ETP` ~ gene)
    )

    saveRDS(result, file=args$outfile)
})
