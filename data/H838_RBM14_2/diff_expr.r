library(dplyr)
library(DESeq2)
library(patchwork)
sys = import('sys')
gset = import('genesets')
plt = import('plot')

de_time = function(.time) {
    eset2 = eset[,eset$time == .time]
    design(eset2) = ~ treatment
    DESeq2::DESeq(eset2) %>%
        DESeq2::results(name = "treatment_RBM14_vs_luc") %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_name") %>%
        as_tibble() %>%
        arrange(padj, pvalue)
}

plot_row = function(..., time) {
    volcs = lapply(list(...), plt$volcano)
    plt$text(time) / wrap_plots(volcs, nrow=1) + plot_layout(heights=c(1,15))
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../../config.yaml'),
        opt('e', 'eset', 'rds', 'eset.rds'),
        opt('o', 'outfile', 'rds', 'diff_expr.rds'),
        opt('p', 'plotfile', 'pdf', 'diff_expr.pdf')
    )

    cfg = yaml::read_yaml(args$config)
    eset = readRDS(args$eset)

    sets = gset$get_human(c("MSigDB_Hallmark_2020", "DoRothEA", "GO_Biological_Process_2021"))

    res = tibble(time = unique(eset$time)) %>%
        rowwise() %>%
        mutate(genes = list(de_time(time)),
#               gsets = do.call(tibble, lapply(sets, function(s) list(gset$test_lm(genes, s))))) %>%
               MSigDB_Hallmark_2020 = list(gset$test_lm(genes, sets$MSigDB_Hallmark_2020)),
               DoRothEA = list(gset$test_lm(genes, sets$DoRothEA)),
               GO_Biological_Process_2021 = list(gset$test_lm(genes, sets$GO_Biological_Process_2021))) %>%
        ungroup() #%>%
#        tidyr::unnest_wider(gsets)

    saveRDS(res, file=args$outfile)

    plots = purrr::pmap(res, plot_row)
    pdf(args$plotfile, 24, 8)
    for (p in plots)
        print(p)
    dev.off()
})
