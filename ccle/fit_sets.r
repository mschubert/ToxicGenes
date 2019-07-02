library(dplyr)
sys = import('sys')
gset = import('data/genesets')

fit_set = function(set, sets, emat, copies, covar=1) {
    genes = intersect(sets[[set]], rownames(emat))
    df = data.frame(expr = c(emat[genes,]),
                    copies = c(copies[genes,]),
                    covar = rep(covar, length(genes))) %>%
        na.omit()
    if (length(unique(na.omit(covar))) > 1)
        fml = expr ~ covar + copies
    else
        fml = expr ~ copies

    mobj = MASS::rlm(fml, data=df, maxit=100)
    mod = broom::tidy(mobj) %>%
        filter(term == "copies") %>%
        select(-term) %>%
        mutate(size = length(genes),
               p.value = sfsmisc::f.robftest(mobj, var="copies")$p.value)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('c', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '6144'),
        opt('o', 'outfile', 'xlsx', 'pan.xlsx'))

    dset = readRDS(args$infile)
    if (args$tissue != "pan")
        dset$idx$tcga_code[dset$idx$tcga_code != args$tissue] = NA
    emat = dset$eset / rowMeans(dset$eset, na.rm=TRUE) - 1

    sets = readRDS(args$setfile) %>%
        gset$filter(min=4, valid=rownames(dset$eset))

    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
        dev = function(x) abs(x-2),
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        tibble(set = names(sets)) %>%
            mutate(res = clustermq::Q(fit_set, set=set, pkgs="dplyr",
                const = list(sets=sets, emat=emat, copies=ff(dset$copies),
                             covar=dset$idx$tcga_code),
                n_jobs=as.integer(args$cores), memory=as.integer(args$memory))) %>%
            tidyr::unnest() %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
})
