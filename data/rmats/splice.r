library(dplyr)
library(patchwork)
library(ggplot2)
plt = import('plot')
gset = import('genesets')
idmap = import('process/idmap')
sys = import('sys')

read_one = function(fname) {
    readr::read_tsv(fname) %>%
        select(gene=geneSymbol, chr, strand, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference) %>%
        mutate(stat = -log10(PValue) * sign(IncLevelDifference)) %>%
        arrange(FDR, PValue)
}

plot_one = function(df, title) {
    plt$volcano(df %>% mutate(circle = gene %in% U12, size=3),
                label="gene", x="IncLevelDifference", y="FDR") +
        ggtitle(title)
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'comp', 'comparison', 'splice-HCC70_rbm8_vs_luc8'),
        opt('j', 'junction', 'JC|JCEC', 'JC'),
        opt('p', 'plotfile', 'pdf', 'splice-rbm8_vs_luc8-JC.pdf')
    )

    U12 = read.table("https://introndb.lerner.ccf.org/static/bed/GRCh38_U12.bed") %>%
        pull(V4) %>%
        sub("GRCh38-([A-Z0-9]+)@.*", "\\1", .) %>%
        idmap$gene(to="hgnc_symbol") %>% unname()
    #U2 = read.table("https://introndb.lerner.ccf.org/static/bed/GRCh38_U2.bed")

    sets = gset$get_human(c("MSigDB_Hallmark_2020", "DoRothEA", "GO_Biological_Process_2021"))

    stypes = c("A3SS", "A5SS", "MXE", "RI", "SE")
    fnames = file.path("rmats_out", args$comp, sprintf("%s.MATS.%s.txt", stypes, args$junction))
    res = tibble(stype=stypes, genes=lapply(fnames, read_one)) %>%
        rowwise() %>%
    #    mutate(GO_Biological_Process_2021 = list(gset$test_lm(genes, sets$GO_Biological_Process_2021))) %>%
        ungroup()

    p1 = mapply(plot_one, res$genes, res$stype, SIMPLIFY=FALSE)
    #p2 = mapply(function(x, n) plt$volcano(x, text.size=2.5) + ggtitle(n),
    #            res$GO_Biological_Process_2021, res$stype, SIMPLIFY=FALSE)

    pdf(args$plotfile, 16, 10)
    print(plt$text("genes") / wrap_plots(p1) + plot_layout(heights=c(1,15)))
    #plt$text("GO_Biological_Process_2021") / wrap_plots(p2) + plot_layout(heights=c(1,15))
    dev.off()
})
