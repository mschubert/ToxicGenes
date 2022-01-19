library(dplyr)
library(patchwork)
library(ggplot2)
plt = import('plot')
idmap = import('process/idmap')
gset = import('genesets')
sys = import('sys')

#' Read differential splicing results from rMATS
#'
#' @param comp      Comparison (common part of file name)
#' @param stypes    Splicing event type (A3SS, MXE, etc.)
#' @param junction  Description of the quantification type (JC, JCEC)
#' @return          A nested tibble containing different results
read_all = function(comp, stypes=c("A3SS", "A5SS", "MXE", "RI", "SE"), junction="JC") {
    read_one = function(fname) {
        readr::read_tsv(fname) %>%
            select(label=GeneID, chr, strand, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference) %>%
            mutate(stat = -log10(PValue) * sign(IncLevelDifference)) %>%
            arrange(FDR, PValue)
    }
    fnames = file.path("rmats_out", comp, sprintf("%s.MATS.%s.txt", stypes, junction))
    tibble(stype=stypes, genes=lapply(fnames, read_one))
}

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
                    x=c("IncLevelDifference", "estimate"), y=c("FDR", "adj.p")) +
            ggtitle(title)
    }
    plots = mapply(plot_one, df$diff_splice, df$stype, SIMPLIFY=FALSE)
    plt$text(sprintf("%s :: %s (%s)", comp, coll, junc)) / wrap_plots(plots) + plot_layout(heights=c(1,15))
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'comp', 'comparison', 'splice-HCC70_rbm8_vs_luc8'),
        opt('o', 'outfile', 'rds', 'splice-rbm8_vs_luc8-JC.rds'),
        opt('p', 'plotfile', 'pdf', 'splice-rbm8_vs_luc8-JC.pdf')
    )

    U12 = read.table("GRCh38_U12.bed") %>%
        pull(V4) %>%
        sub("GRCh38-([A-Z0-9]+)@.*", "\\1", .) %>%
        idmap$gene(to="hgnc_symbol") %>% unname()
    #U2 = read.table("GRCh38_U2.bed")

    sets = gset$get_human(c("MSigDB_Hallmark_2020", "DoRothEA", "GO_Biological_Process_2021"))

    res = list(
        jc = read_all(args$comp, junction="JC"),
        jcec = read_all(args$comp, junction="JCEC")
    ) %>% bind_rows(.id="junction") %>%
        rowwise() %>%
            mutate(MSigDB_Hallmark_2020 = list(gset$test_lm(genes, sets$MSigDB_Hallmark_2020)),
                   DoRothEA = list(gset$test_lm(genes, sets$DoRothEA)),
                   GO_Biological_Process_2021 = list(gset$test_lm(genes, sets$GO_Biological_Process_2021))) %>%
        ungroup()

    plots = res %>%
        tidyr::gather("collection", "diff_splice", -junction, -stype) %>%
        group_by(junction, collection) %>%
            tidyr::nest() %>%
        rowwise() %>%
            mutate(plot = list(plot_asm(args$comp, collection, junction, data, U12))) %>%
        ungroup()

    pdf(args$plotfile, 16, 10)
    for (p in plots$plot)
        print(p)
    dev.off()

    saveRDS(res, file=args$outfile)
})
