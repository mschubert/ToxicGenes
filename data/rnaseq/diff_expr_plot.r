library(dplyr)
library(patchwork)
library(ggplot2)
plt = import('plot')
idmap = import('process/idmap')
sys = import('sys')

#' Assemble panels for different splice events
#'
#' @param coll  Name of collection (genes, MSigDB_Hallmark, etc.)
#' @param df    A data.frame with fields: stype [A3SS, MXE, etc.], diff_splice [df]
#' @param hl    Character vector of genes to highlight with circles
#' @return      A patchwork object of volcano plots
plot_asm = function(coll, df, hl=c()) {
    plot_one = function(df, title) {
        plt$volcano(df, size = c("size_used", "baseMean")) +
            ggtitle(title)
    }
    plots = mapply(plot_one, df$diff_expr, df$cond, SIMPLIFY=FALSE)
    plt$text(coll) / wrap_plots(plots, byrow=FALSE) + plot_layout(heights=c(1,15))
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'diff_expr.rds'),
        opt('p', 'plotfile', 'pdf', 'diff_expr.pdf')
    )

    res = readRDS(args$infile)

    plots = res %>%
        tidyr::gather("collection", "diff_expr", -cond) %>%
        group_by(collection) %>%
            tidyr::nest() %>%
        rowwise() %>%
            mutate(plot = list(plot_asm(collection, data))) %>%
        ungroup()

    pdf(args$plotfile, 16, 10)
    for (p in plots$plot)
        print(p)
    dev.off()
})
