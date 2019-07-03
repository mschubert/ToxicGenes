library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fit = function(genes, emat, copies, covar=1) {
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
        mutate(n_aneup = sum(abs(df$cancer_copies-2) > 0.2),
               n_genes = length(genes),
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

    #TODO: check where these errors come from
    if (args$tissue == "SKCM") {
        emat = emat[rownames(emat) != "LINC00320",]
        dset$copies = dset$copies[rownames(dset$copies) != "LINC00320",]
    }
    if (args$tissue == "BRCA") {
        emat = emat[! rownames(emat) %in% c("DNAJA1P5","DEFA4"),]
        dset$copies = dset$copies[! rownames(dset$copies) %in% c("DNAJA1P5","DEFA4"),]
    }

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
                             covar=dset$idx$tcga_code)) %>%
            setNames(names(sets)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
})
