library(dplyr)
library(patchwork)
library(ggplot2)
plt = import('plot')
idmap = import('process/idmap')
sys = import('sys')

read_one = function(fname) {
    readr::read_tsv(fname) %>%
        select(gene=geneSymbol, chr, strand, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference) %>%
        mutate(stat = -log10(PValue) * sign(IncLevelDifference)) %>%
        arrange(FDR, PValue)
}

read_all = function(comp, stypes=c("A3SS", "A5SS", "MXE", "RI", "SE"), junction="JC") {
    fnames = file.path("rmats_out", comp, sprintf("%s.MATS.%s.txt", stypes, junction))
    res = tibble(stype=stypes, genes=lapply(fnames, read_one))
}

plot_one = function(df, title) {
    plt$volcano(df %>% mutate(circle = gene %in% U12, size=3),
                label="gene", x="IncLevelDifference", y="FDR") +
        ggtitle(title)
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'comp', 'comparison', 'splice-HCC70_rbm8_vs_luc8'),
        opt('p', 'plotfile', 'pdf', 'splice-rbm8_vs_luc8-JC.pdf')
    )

    U12 = read.table("GRCh38_U12.bed") %>%
        pull(V4) %>%
        sub("GRCh38-([A-Z0-9]+)@.*", "\\1", .) %>%
        idmap$gene(to="hgnc_symbol") %>% unname()
    #U2 = read.table("GRCh38_U2.bed")

    jc = read_all(args$comp, junction="JC")
    jcec = read_all(args$comp, junction="JCEC")

    p1 = mapply(plot_one, jc$genes, jc$stype, SIMPLIFY=FALSE)
    p2 = mapply(plot_one, jcec$genes, jcec$stype, SIMPLIFY=FALSE)

    pdf(args$plotfile, 16, 10)
    print(plt$text(sprintf("%s (JC)", args$comp)) / wrap_plots(p1) + plot_layout(heights=c(1,15)))
    print(plt$text(sprintf("%s (JCEC)", args$comp)) / wrap_plots(p2) + plot_layout(heights=c(1,15)))
    dev.off()
})
