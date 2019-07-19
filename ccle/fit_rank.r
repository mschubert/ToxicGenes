library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fit = function(genes, emat, copies, covar=1) {
    cov_noNA = covar
    cov_noNA[is.na(covar)] = "unknown"
    crank = narray::map(copies[genes,,drop=FALSE], along=2, subsets=cov_noNA,
                        function(x) rank(x) / length(x) - 0.5)
    erank = narray::map(emat[genes,,drop=FALSE], along=2, subsets=cov_noNA,
                        function(x) rank(x) / length(x) - 0.5)
    df = data.frame(expr = c(emat[genes,,drop=FALSE] /
                             rowMeans(emat[genes,,drop=FALSE], na.rm=TRUE)),
                    copies = c(copies[genes,]),
                    erank = c(erank),
                    crank = c(crank),
                    covar = rep(covar, length(genes))) %>%
        na.omit() %>%
        sample_n(min(nrow(.), 1e5))
    if (length(unique(na.omit(covar))) > 1)
        fml = erank ~ covar + crank
    else
        fml = erank ~ crank

    mod = lm(fml, data=df) %>%
        broom::tidy() %>%
        filter(term == "crank") %>%
        select(-term) %>%
        mutate(n_aneup = sum(abs(df$copies-2) > 0.2),
               n_genes = length(genes))
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('c', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '6144'),
        opt('o', 'outfile', 'xlsx', 'pan_rank/genes.xlsx'))

    dset = readRDS(args$infile)
    if (args$tissue != "pan")
        dset$clines$tcga_code[dset$clines$tcga_code != args$tissue] = NA

    emat = dset$eset # already copy-normalized in dset
    if (grepl("genes\\.xlsx", args$outfile))
        sets = setNames(rownames(emat), rownames(emat))
    else
        sets = readRDS(args$setfile) %>%
            gset$filter(min=4, valid=rownames(emat))

    w = clustermq::workers(n_jobs = as.integer(args$cores),
                           template = list(memory = as.integer(args$memory)))

    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
        dev = function(x) abs(x-2),
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        res = clustermq::Q(do_fit, genes=sets, workers=w, pkgs="dplyr",
                const = list(emat=emat, copies=ff(dset$copies),
                             covar=dset$clines$tcga_code)) %>%
            setNames(names(sets)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
    w$cleanup()
})
