library(dplyr)
library(ggplot2)
plt = import('plot')
sys = import('sys')
gset = import('data/genesets')

do_fit = function(sets, expr) {
    expr %>%
        mutate(set = `GENE SYMBOL` %in% sets) %>%
        lm(z_LFC ~ set, data=.) %>%
        broom::tidy() %>%
        filter(term == "setTRUE") %>%
        select(-term) %>%
        mutate(size = length(sets))
}

do_plot = function(res) {
    res$label = res$name
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
        opt('t', 'tissue', 'pan|TCGA', 'pan'),
        opt('o', 'outfile', 'xls', 'pan/CH.HALLMARK.xlsx'),
        opt('p', 'plotfile', 'pdf', 'pan/CH.HALLMARK.pdf'))

    expr = readRDS(args$infile)

    if (grepl("genes\\.xlsx", args$outfile)) {
        genes = unique(expr$`GENE SYMBOL`)
        sets = setNames(as.list(genes), genes)
    } else
        sets = readRDS(args$setfile) %>%
            gset$filter(min=4, valid=expr$`GENE SYMBOL`)

    if (args$tissue != "pan")
        expr = filter(expr, tissue %in% args$tissue)

    result = clustermq::Q(do_fit, sets=sets, n_jobs=5, pkgs="dplyr",
            const=list(expr=expr)) %>%
        setNames(names(sets)) %>%
        bind_rows(.id="name") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)

    title = paste(args$tissue, tools::file_path_sans_ext(basename(args$setfile)))
    pdf(args$plotfile)
    print(do_plot(result) + ggtitle(title))
    dev.off()

    writexl::write_xlsx(result, args$outfile)
})
