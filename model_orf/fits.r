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
        opt('f', 'field', 'response variable', 'LFC DMSO/ETP'),
        opt('o', 'outfile', 'xls', 'fits_naive.xlsx'),
        opt('p', 'plotfile', 'pdf', 'fits_naive.pdf'))

    expr = readRDS(args$infile) %>%
        mutate(`LFC DMSO/ETP` = `LFC DMSO/ETP` + runif(nrow(.)) * 0.01)

    pan = list(
        pan = do_fit(expr, as.formula(sprintf("`%s` ~ gene", args$field))),
        pancov = do_fit(expr, as.formula(sprintf("`%s` ~ tissue + gene", args$field)))
    )

    tissue = sapply(c("NB", "OV", "BRCA", "SKCM", "MB"),
                    function(t) expr %>%
                        filter(tissue %in% t) %>%
                        do_fit(as.formula(sprintf("`%s` ~ gene", args$field))),
                    simplify=FALSE)

    result = c(pan, tissue)

    pdf(args$plotfile)
    for (i in seq_along(result)) {
        cur = result[[i]] %>% dplyr::rename(label = `GENE SYMBOL`)
        plt$volcano(cur, label_top=30, repel=TRUE, x_label_bias=4, p=0.1)
    }
    dev.off()

    writexl::write_xlsx(result, args$outfile)
})
