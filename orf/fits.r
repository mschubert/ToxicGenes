library(dplyr)
library(cowplot)
plt = import('plot')
sys = import('sys')

do_fit = function(expr, fml) {
    ffun = function(g) {
        expr %>%
            mutate(gene = `Construct IDs` == g) %>%
            lm(fml, data=.) %>%
            broom::tidy()
    }

    res = expr %>%
        select(`GENE SYMBOL`, `Construct IDs`) %>%
        distinct() %>%
        mutate(result = clustermq::Q(ffun, `Construct IDs`, n_jobs=10,
            export=list(expr=expr, fml=fml), pkgs="dplyr")) %>%
        tidyr::unnest() %>%
        filter(term == "geneTRUE") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

do_plot = function(res) {
    res$label = res$`GENE SYMBOL`
    if (all(res$estimate[rank(res$p.value) < 10] > 0))
        res$label[res$estimate > 0] = NA

    res %>%
        mutate(size=5) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", thresh=0.1, dir=-1) %>%
        plt$volcano(label_top=30, repel=TRUE, x_label_bias=4)
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
        pan = do_fit(expr, `LFC DMSO/ETP` ~ gene),
        pancov = do_fit(expr, `LFC DMSO/ETP` ~ tissue + gene)
    )

    tissue = sapply(c("NB", "OV", "BRCA", "SKCM", "MB"),
                    function(t) expr %>%
                        filter(tissue %in% t) %>%
                        do_fit(`LFC DMSO/ETP` ~ gene),
                    simplify=FALSE)

    result = c(pan, tissue)

    pdf(args$plotfile)
    for (i in seq_along(result))
        print(do_plot(result[[i]]) + ggtitle(names(result)[i]))
    dev.off()

    writexl::write_xlsx(result, args$outfile)
})
