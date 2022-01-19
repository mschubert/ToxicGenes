library(dplyr)
library(patchwork)
library(ggplot2)
plt = import('plot')
idmap = import('process/idmap')
sys = import('sys')

#' Assemble panels for different splice events
#'
#' @param comp  Name of comparison (common label for all plots)
#' @param coll  Name of collection (genes, MSigDB_Hallmark, etc.)
#' @param junc  Splice quantification type (JC, JCEC)
#' @param df    A data.frame with fields: stype [A3SS, MXE, etc.], diff_splice [df]
#' @param hl    Character vector of genes to highlight with circles
#' @return      A patchwork object of volcano plots
plot_asm = function(comp, coll, junc, df, hl=c()) {
    plot_one = function(df, title) {
        plt$volcano(df %>% mutate(circle = label %in% hl, size=3),
                    size = c("size_used", "size"),
                    x=c("IncLevelDifference", "estimate"), y=c("FDR", "adj.p")) +
            ggtitle(title)
    }
    plots = mapply(plot_one, df$diff_splice, df$stype, SIMPLIFY=FALSE)
    plt$text(sprintf("%s :: %s (%s)", comp, coll, junc)) / wrap_plots(plots) + plot_layout(heights=c(1,15))
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'splice/splice-H1650_rbm8_vs_luc8.rds'),
        opt('p', 'plotfile', 'pdf', 'splice/splice-H1650_rbm8_vs_luc8.pdf')
    )

    res = readRDS(args$infile)
    comp = tools::file_path_sans_ext(basename(args$infile))

    U12 = read.table("GRCh38_U12.bed") %>%
        pull(V4) %>%
        sub("GRCh38-([A-Z0-9]+)@.*", "\\1", .) %>%
        idmap$gene(to="hgnc_symbol") %>% unname()
    #U2 = read.table("GRCh38_U2.bed")

    plots = res %>%
        tidyr::gather("collection", "diff_splice", -junction, -stype) %>%
        group_by(junction, collection) %>%
            tidyr::nest() %>%
        rowwise() %>%
            mutate(plot = list(plot_asm(comp, collection, junction, data, U12))) %>%
        ungroup()

    pdf(args$plotfile, 16, 10)
    for (p in plots$plot)
        print(p)
    dev.off()
})
