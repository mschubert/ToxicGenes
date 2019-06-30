library(dplyr)
library(cowplot)
plt = import('plot')
sys = import('sys')

do_fit = function(expr, sets) {
    ffun = function(s) {
        gset = sets[[s]]
        expr %>%
            mutate(set = `GENE SYMBOL` %in% gset) %>%
            lm(z_LFC ~ set, data=.) %>%
            broom::tidy() %>%
            mutate(size = length(gset))
    }

    res = tibble(set = names(sets)) %>%
        mutate(result = clustermq::Q(ffun, set, n_jobs=10,
            export=list(expr=expr, sets=sets), pkgs="dplyr")) %>%
        tidyr::unnest() %>%
        filter(term == "setTRUE") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

do_plot = function(res) {
    res$label = res$set
    if (all(res$estimate[rank(res$p.value) < 10] > 0))
        res$label[res$estimate > 0] = NA

    res %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate", thresh=0.1, dir=-1) %>%
        plt$volcano(base.size=0.2, label_top=30, repel=TRUE, x_label_bias=4)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'overview.rds'),
        opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
        opt('o', 'outfile', 'xls', 'sets/CH.HALLMARK.xlsx'),
        opt('p', 'plotfile', 'pdf', 'sets/CH.HALLMARK.pdf'))

    expr = readRDS(args$infile)
    sets = readRDS(args$setfile)

    pan = do_fit(expr, sets)
    tissue = sapply(c("NB", "OV", "BRCA", "SKCM", "MB"),
                    function(t) do_fit(filter(expr, tissue %in% t), sets),
                    simplify=FALSE)
    result = c(list(pan=pan), tissue)

    pdf(args$plotfile)
    for (i in seq_along(result))
        print(do_plot(result[[i]]) + ggtitle(names(result)[i]))
    dev.off()

    writexl::write_xlsx(result, args$outfile)
})
