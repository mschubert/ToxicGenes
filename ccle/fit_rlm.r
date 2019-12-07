library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fit = function(genes, emat, copies, covar=1, et=0.15) {
    df = data.frame(expr = c(emat[genes,,drop=FALSE] /
                             rowMeans(emat[genes,,drop=FALSE], na.rm=TRUE)),
                    copies = c(copies[genes,]),
                    covar = rep(covar, length(genes))) %>%
        na.omit() %>%
        sample_n(min(nrow(.), 1e5))
    if (length(unique(na.omit(covar))) > 1)
        fml = expr ~ covar + copies
    else
        fml = expr ~ copies

    tryCatch({
        mobj = MASS::rlm(fml, data=df, maxit=100)
        mod = broom::tidy(mobj) %>%
            filter(term == "copies") %>%
            select(-term) %>%
            mutate(n_aneup = sum(abs(df$copies-2) > et),
                   n_genes = length(genes),
                   p.value = sfsmisc::f.robftest(mobj, var="copies")$p.value)
    }, error = function(e) {
        warning(genes, " : ", conditionMessage(e), immediate.=TRUE)
        data.frame(estimate=NA)
    })
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('c', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '6144'),
        opt('o', 'outfile', 'xlsx', 'pan_rlm/genes.xlsx'))

    et = yaml::read_yaml(args$config)$euploid_tol

    dset = readRDS(args$infile)
    if (args$tissue != "pan")
        dset$clines$tcga_code[dset$clines$tcga_code != args$tissue] = NA

    emat = dset$eset # already copy-normalized in dset
    genes = setNames(rownames(emat), rownames(emat))

    w = clustermq::workers(n_jobs = as.integer(args$cores),
                           template = list(memory = as.integer(args$memory)))

    ffuns = list(
        amp = function(x) { x[x < 2-et] = NA; x },
        del = function(x) { x[x > 2+et] = NA; x },
        dev = function(x) abs(x-2),
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        res = clustermq::Q(do_fit, genes=genes, workers=w, pkgs="dplyr",
                const = list(emat=emat, copies=ff(dset$copies),
                             covar=dset$clines$tcga_code, et=et)) %>%
            setNames(names(genes)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
    w$cleanup()
})
