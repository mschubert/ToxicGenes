library(dplyr)
sys = import('sys')
plt = import('plot')
gset = import('data/genesets')

do_fit = function(set, sets, emat, copies, covar=1) {
    genes = intersect(sets[[set]], rownames(emat))
    df = data.frame(expr = c(emat[genes,] / rowMeans(emat[genes,], na.rm=TRUE) - 1),
                    copies = c(copies[genes,]),
                    covar = rep(covar, length(genes))) %>%
        na.omit()

    has_covar = length(unique(na.omit(covar))) > 1
    if (has_covar) {
        fml = expr ~ covar + copies
    } else {
        fml = expr ~ copies
    }

    mobj = MASS::rlm(fml, data=df, maxit=100)
    mod = broom::tidy(mobj) %>%
        filter(term == "copies") %>%
        select(-term) %>%
        mutate(size = length(genes),
               p.value = sfsmisc::f.robftest(mobj, var="copies")$p.value)
}

all_fits = function(sets, emat, copies, tissues=NA) {
    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
        dev = function(x) abs(x-2),
        all = identity
    )

    do_ffun = function(ffun) {
        tibble(set = names(sets)) %>%
            mutate(res = clustermq::Q(do_fit, set=set, n_jobs=3,
                const = list(sets=sets, emat=emat, copies=ffun(copies), covar=tissues),
                pkgs = "dplyr")) %>%
            tidyr::unnest() %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    }
    lapply(ffuns, do_ffun)
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
    sets = readRDS(args$setfile) %>%
        gset$filter(min=4, valid=rownames(dset$eset))

    fits = all_fits(sets, dset$eset, dset$copies, dset$idx$tcga_code)
    plots = lapply(fits, do_plot)

    pdf(args$plotfile)
    for (i in seq_along(plots))
        print(plots[[i]] + ggtitle(names(plots)[i]))
    dev.off()

    writexl::write_xlsx(fits, args$outfile)
})
