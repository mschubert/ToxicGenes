library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fit = function(genes, emat, copies, covar=1, et=0.15) {
    df = data.frame(gene = rep(genes, each=ncol(emat)),
                    expr = c(emat[genes,,drop=FALSE]),
                    copies = c(copies[genes,]),
                    covar = rep(covar, length(genes))) %>%
#        filter(copies >= 0.5) %>% # otherwise genes with 0 copies and FP reads
        mutate(expr = expr * copies / 2) %>% #, # undo normmatrix normalization
#        filter(expr < quantile(expr, 0.95) & copies < quantile(copies, 0.99)) %>%
        na.omit() %>%
        sample_n(min(nrow(.), 1e5))
    if (length(unique(na.omit(covar))) > 1)
        fml = expr - copies ~ covar + copies
    else
        fml = expr - copies ~ copies

    tryCatch({
        # fit intercept model
        df2 = df %>%
            filter(copies > 2-et & copies < 2+et) %>%
            group_by(gene, covar) %>%
            summarize(intcp = MASS::rlm(expr ~ 1, maxit=100)$coefficients["(Intercept)"]) %>%
            inner_join(df, by=c("gene", "covar")) %>%
            mutate(expr = expr / intcp - 1,
                   copies = copies / 2 - 1) %>%
            na.omit()

        # fit actual model
        mobj = MASS::rlm(fml, data=df2, maxit=100)
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
        opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('j', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '6144'),
        opt('o', 'outfile', 'xlsx', 'pan_rlm/genes.xlsx'))

    et = yaml::read_yaml(args$config)$euploid_tol

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
        amp = function(x) { x[x < 2-et] = NA; x },
        del = function(x) { x[x > 2+et] = NA; x },
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        res = clustermq::Q(do_fit, genes=sets, workers=w, pkgs="dplyr",
                const = list(emat=emat, copies=ff(dset$copies),
                             covar=dset$clines$tcga_code, et=et)) %>%
            setNames(names(sets)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
    w$cleanup()
})
