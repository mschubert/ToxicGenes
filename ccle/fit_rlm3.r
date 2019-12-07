library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fit = function(genes, emat, copies, covar=1, et=0.15) {
    df = data.frame(gene = rep(genes, each=ncol(emat)),
                    expr = c(emat[genes,,drop=FALSE]),
                    copies = c(copies[genes,]),
                    covar = rep(covar, length(genes))) %>%
#        filter(copies >= 0.5) %>% # otherwise genes with 0 copies and FP reads
        mutate(covar = paste(covar, gene, sep=":"),
               expr = expr * copies / 2) %>% #, # undo normmatrix normalization
#        filter(expr < quantile(expr, 0.95) & copies < quantile(copies, 0.99)) %>%
        na.omit() %>%
        sample_n(min(nrow(.), 1e5))

    tryCatch({
        if (length(unique(na.omit(covar))) > 1) {
            fml = expr ~ covar + copies
            fml_mean = expr ~ covar
            mobj = MASS::rlm(fml, data=df, maxit=100)
            ucovar = unique(df$covar)
            pred = predict(mobj, newdata=expand.grid(covar=ucovar, copies=2))
            expr_per_copy = data.frame(covar=ucovar, expr_per_copy=pred/2)
            df = inner_join(df, expr_per_copy, by="covar")
        } else {
            fml = expr ~ copies
            fml_mean = expr ~ 1
            mobj = MASS::rlm(fml, data=df, maxit=100)
            pred = predict(mobj, newdata=data.frame(copies=2))
            df$expr_per_copy = pred / 2
        }

        if (mean(pred) < 0) # maybe add check for mobj2/intercept as well?
            stop("predicted negative expression for euploid")

        df$expr = with(df, expr - copies * expr_per_copy)
        mobj2 = MASS::rlm(fml, data=df, maxit=100)
        mmean = MASS::rlm(fml_mean, data=df, maxit=100)
        res = broom::tidy(mobj2) %>%
            filter(term == "copies") %>%
            select(-term) %>%
            mutate(estimate = 2 * estimate / mean(pred), # pct_comp
                   n_aneup = sum(abs(df$copies-2) > et),
                   n_genes = length(genes),
                   eup_reads = mean(pred),
                   slope_diff = 2 * estimate,
                   rsq = 1 - sum(mobj2$w * mobj2$resid^2) / sum(mmean$w * mmean$resid^2),
                   p.value = sfsmisc::f.robftest(mobj2, var="copies")$p.value)

    }, error = function(e) {
        warning(genes, ": ", conditionMessage(e), immediate.=TRUE)
        data.frame(estimate=NA)
    })
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('j', 'cores', 'integer', '10'),
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
