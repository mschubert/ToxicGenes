library(dplyr)
library(DESeq2)
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')
gset = import('genesets')

process = function(deseq, extract) {
    DESeq2::results(deseq, name=extract) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_name") %>%
        as_tibble() %>%
        arrange(padj, pvalue) %>%
        select(gene_name, everything())
}

de_all = function(eset) {
    design(eset) = ~ cond
    res = DESeq2::DESeq(eset)
    process(res, "cond_RBM14_vs_Luc")
}

de_all_covar = function(eset) {
    design(eset) = ~ cline + cond
    res = DESeq2::DESeq(eset)
    process(res, "cond_RBM14_vs_Luc")
}

de_cline = function(eset, cline) {
    eset2 = eset[,eset$cline == cline]
    de_all(eset2)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'eset.rds'),
        opt('t', 'time', '8h|24h|all', '8h'),
        opt('o', 'outfile', 'rds', 'diff_expr.rds')
    )

    eset = readRDS(args$infile)
    eset = eset[,eset$time != "72h"] # only H838
    if (args$time != "all")
        eset = eset[,eset$time == args$time]

    sets = gset$get_human(c("MSigDB_Hallmark_2020", "CORUM_all", "CORUM_core", "CORUM_splice",
                            "HMDB_Metabolites", "miRTarBase_2017", "DoRothEA", "GO_Biological_Process_2021",
                            "GO_Cellular_Component_2021", "GO_Molecular_Function_2021"), conf="A") %>%
        gset$filter(max=500)

    de_genes = c(list(all=de_all(eset), all_covar=de_all_covar(eset)),
                 sapply(levels(eset$cline), de_cline, eset=eset, simplify=FALSE))

    res = tibble(
        cond = names(de_genes),
        genes = unname(de_genes)
    ) %>% rowwise()

    for (sname in names(sets))
        res = res %>%
            mutate({{ sname }} := list(gset$test_lm(genes, sets[[sname]])))

    saveRDS(ungroup(res), file=args$outfile)
})
