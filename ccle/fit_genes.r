library(dplyr)
sys = import('sys')

fit_gene = function(expr, copies, covar=1, ffun=identity) {
    df = na.omit(data.frame(expr=expr, copies=ffun(copies), covar=covar))
    if (length(unique(na.omit(covar))) > 1)
        fml = expr ~ covar + copies
    else
        fml = expr ~ copies

    mobj = MASS::rlm(fml, data=df, maxit=100)
    mod = broom::tidy(mobj) %>%
        filter(term == "copies") %>%
        select(-term) %>%
        mutate(n_aneup = sum(abs(df$copies-2) > 0.2),
               size = n_aneup,
               p.value = sfsmisc::f.robftest(mobj, var="copies")$p.value)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('c', 'cores', 'integer', '10'), # defunct
        opt('m', 'memory', 'integer', '6144'), # defunct
        opt('o', 'outfile', 'xlsx', 'pan.xlsx'))

    #clustermq::register_dopar_cmq(n_jobs=as.integer(args$cores), memory=as.integer(args$memory))

    dset = readRDS(args$infile)
    if (args$tissue != "pan")
        dset$idx$tcga_code[dset$idx$tcga_code != args$tissue] = NA
    emat = dset$eset / rowMeans(dset$eset, na.rm=TRUE) - 1
    cmat = dset$copies

    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
        dev = function(x) abs(x-2),
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        narray::lambda(~ fit_gene(emat, cmat, dset$idx$tcga_code, ff),
                       along = c(emat=1, cmat=1),
                       simplify=FALSE, expand_grid=FALSE) %>%
            tidyr::unnest() %>%
            rename(gene = emat) %>%
            select(-cmat) %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
})
