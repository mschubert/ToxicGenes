library(dplyr)
library(ggplot2)
library(patchwork)
plt = import('plot')
idmap = import('process/idmap')
sys = import('sys')

#' Assemble panels for different splice events
#'
#' @param junction  Splice quantification type (JC, JCEC)
#' @param stype     Splice ftype (A3SS, A5SS, SE, RI, MXE)
#' @param df        A data.frame of gene set gene set differences
#' @return          A volcano plot with title
plot_one = function(junction, stype, df, hl=c()) {
    if (nrow(df) == 0)
        return(plt$text("No observations"))
    plt$volcano(df %>% mutate(circle = label %in% hl), size = c("size_used", "size"),
                x=c("IncLevelDifference", "mean_shrunkAbsIncDiff"), y=c("FDR", "adj.p")) +
        ggtitle(sprintf("%s (%s)", stype, junction))
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
        tidyr::pivot_longer(-c(junction, stype)) %>%
        mutate(name = factor(name, levels=unique(name))) %>%
        rowwise() %>%
            mutate(plot = list(plot_one(junction, stype, value, U12))) %>%
        group_by(name) %>%
            summarize(asm = list((plt$text(name[1], size=8) / wrap_plots(plot, nrow=2)) +
                                 plot_layout(heights=c(1,20))))

    pdf(args$plotfile, 22, 10)
    for (p in plots$asm)
        print(p)
    dev.off()
})
