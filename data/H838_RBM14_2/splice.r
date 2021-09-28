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
        arrange(FDR, PValue)
}

plot_one = function(df, title) {
    plt$volcano(df %>% mutate(circle = gene %in% U12, size=3),
                label="gene", x="IncLevelDifference", y="FDR") +
        ggtitle(title)
}

args = sys$cmd$parse(
    opt('p', 'plotfile', 'pdf', 'splice.pdf')
)

U12 = read.table("https://introndb.lerner.ccf.org/static/bed/GRCh38_U12.bed") %>%
    pull(V4) %>%
    sub("GRCh38-([A-Z0-9]+)@.*", "\\1", .) %>%
    idmap$gene(to="hgnc_symbol") %>% unname()
#U2 = read.table("https://introndb.lerner.ccf.org/static/bed/GRCh38_U2.bed")

sets = gset$get_human(c("MSigDB_Hallmark_2020", "DoRothEA", "GO_Biological_Process_2021"))

stypes = c("A3SS", "A5SS", "MXE", "RI", "SE")
fnames = file.path("rmats_out/rbm8_vs_luc8", sprintf("%s.MATS.JC.txt", stypes))
res = tibble(stype=stypes, genes=lapply(fnames, read_one)) %>%
    rowwise() %>%
    mutate(GO_Biological_Process_2021 = list(gset$test_lm(genes, sets$GO_Biological_Process_2021,
                                                          stat="IncLevelDifference")))

plots = mapply(plot_one, res$genes, res$stype, SIMPLIFY=FALSE)
pdf(args$plotfile, 16, 10)
wrap_plots(plots)
dev.off()
