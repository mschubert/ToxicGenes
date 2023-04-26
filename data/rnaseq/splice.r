library(dplyr)
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
        message(fname)
        try({
            readr::read_tsv(fname) %>%
                select(label=GeneID, chr, strand, PValue, FDR, IncLevel1, IncLevel2, IncLevelDifference) %>%
                mutate(shrunkIncDiff = (1-(pmin(1,PValue*10))) * IncLevelDifference) %>%
                arrange(FDR, PValue)
        })
    }
    fnames = file.path("rmats_out", comp, sprintf("%s.MATS.%s.txt", stypes, junction))
    re = tibble(stype=stypes, genes=lapply(fnames, read_one))
    re[sapply(re$genes, function(x) class(x)[1]) != "try-error",]
}

do_test = function(genes) {
    tryCatch({
        g = genes %>% group_by(label) %>%
            summarize(shrunkAbsIncDiff = max(abs(shrunkIncDiff)))
        gset$test_lm(g, sets[[sname]], stat="shrunkAbsIncDiff") %>%
            dplyr::rename(mean_shrunkAbsIncDiff=estimate)
    }, error = function(e) tibble())
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'comp', 'comparison', 'splice-HCC70_rbm8_vs_luc8'),
        opt('o', 'outfile', 'rds', 'splice/splice-HCC70_rbm8_vs_luc8.rds')
    )

    sets = gset$get_human(c("MSigDB_Hallmark_2020", "CORUM_all", "CORUM_core", "CORUM_splice",
                            "HMDB_Metabolites", "miRTarBase_2017", "DoRothEA", "GO_Biological_Process_2021",
                            "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"), conf="A") %>%
        gset$filter(max=500)

    res = list(
        jc = read_all(args$comp, junction="JC"),
        jcec = read_all(args$comp, junction="JCEC")
    ) %>% bind_rows(.id="junction")

    for (sname in names(sets))
        res = rowwise(res) %>% mutate({{ sname }} := list(do_test(genes)))

    saveRDS(ungroup(res), file=args$outfile)
})
