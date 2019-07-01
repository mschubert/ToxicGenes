library(dplyr)
sys = import('sys')
plt = import('plot')
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

do_plot = function(data) {
    data %>%
        mutate(label=set) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
        plt$volcano(base.size=0.2, label_top=50, repel=TRUE,
                    text.size=2.5, x_label_bias=5, pos_label_bias=0.15)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('o', 'outfile', 'xlsx', 'pan.xlsx'),
        opt('p', 'plotfile', 'pdf', 'pan.pdf'))

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
            mutate(res = clustermq::Q(fit_set, set=set, n_jobs=3,
                const = list(sets=sets, emat=emat, copies=ff(dset$copies),
                             covar=dset$idx$tcga_code),
                pkgs = "dplyr")) %>%
            tidyr::unnest() %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    pdf(args$plotfile)
    for (i in seq_along(fits))
        print(do_plot(fits[[i]]) + ggtitle(names(fits)[i]))
    dev.off()

    writexl::write_xlsx(fits, args$outfile)
})
