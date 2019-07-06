library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fit = function(genes, emat, copies, covar=1) {
    crank = narray::map(copies[genes,,drop=FALSE], along=2, subsets=covar,
                        function(x) rank(x) / length(x) - 0.5)
    df = data.frame(expr = c(emat[genes,]),
                    copies = c(copies[genes,]),
                    crank = c(crank),
                    covar = rep(covar, length(genes))) %>%
        na.omit() %>%
        sample_n(min(nrow(.), 1e5))
    if (length(unique(na.omit(covar))) > 1)
        fml = expr ~ covar + crank
    else
        fml = expr ~ crank

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
        opt('o', 'outfile', 'xlsx', 'pan.xlsx'))

    dset = readRDS(args$infile)
    if (args$tissue != "pan")
        dset$idx$tcga_code[dset$idx$tcga_code != args$tissue] = NA
#    emat = dset$eset / rowMeans(dset$eset, na.rm=TRUE) - 1

#    copies = narray::like(0, dset$copies)
#    copies[dset$copies < 1.8] = -rank(abs(dset$copies[dset$copies < 1.8]))
#    copies[dset$copies > 2.2] = rank(dset$copies[dset$copies > 2.2])

    # ranks by cohort and gene
    emat = narray::map(dset$eset, along=2, subsets=dset$idx$tcga_code,
                       function(x) rank(x) / length(x) - 0.5)
    dset$copies = dset$copies[,colnames(emat)]
    dset$idx = dset$idx[match(colnames(emat), dset$idx$CCLE_ID),]

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
    w$cleanup()
})
