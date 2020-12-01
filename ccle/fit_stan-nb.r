library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fit = function(genes, eset, copies, covar=1, et=0.15) {
    if (length(unique(covar)) == 1) {
        full = expr ~ eup_equiv
        red = expr ~ 1
    } else {
        full = expr ~ covar + eup_equiv
        red = expr ~ covar
    }

    tryCatch({
        res = rstanarm::stan_glm(full, data=df, offset=sf,
                                 family = neg_binomial_2(link="identity"),
                                 prior = normal(mean(df$expr), sd(df$expr)),
                                 prior_intercept = normal(mean(df$expr), sd(df$expr)),
                                 seed = 2380719)

        rmat = as.matrix(res)
        hdist = statip::hellinger(rmat[,"(Intercept)"], rmat[,"eup_equiv"])
        #for rsq, could compare bayes factors of model +/- eup_equiv

#        res2 = rstanarm::stan_glm(red, data=df, offset=sf,
#                                  family = neg_binomial_2(link="identity"),
#                                  prior = normal(mean(df$expr), sd(df$expr)),
#                                  prior_intercept = normal(mean(df$expr), sd(df$expr)),
#                                  seed = 2380719)
#
#        loo1 = loo(res, k_threshold = 0.7)
#        loo2 = loo(res2, k_threshold = 0.7)
#        cmp = loo_compare(loo1, loo2)

        pseudo_p = pnorm(abs((mean(rmat[,"(Intercept)"]) - mean(rmat[,"eup_equiv"])) /
                             sd(rmat[,"(Intercept)"])), lower.tail=F)

        tibble(estimate = mean(rmat[,"eup_equiv"]) / mean(rmat[,"(Intercept)"]) - 1, # pct_comp
               n_aneup = sum(abs(df$copies-2) > 1-et),
               n_genes = length(genes),
               eup_reads = mean(rmat[,"(Intercept)"]),
               slope_diff = mean(rmat[,"(Intercept)"]) - mean(rmat[,"eup_equiv"]),
               rsq = hdist, # not rsq, but [0,1]
               p.value = pseudo_p)
#               loo_z = -cmp["res2", "elpd_diff"] / cmp["res2", "se_diff"])

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
        opt('o', 'outfile', 'xlsx', 'pan_stan-nb/genes.xlsx'))

    et = yaml::read_yaml(args$config)$euploid_tol

    dset = readRDS(args$infile)
    if (args$tissue != "pan") {
        keep = dset$clines$tcga_code == args$tissue
        dset$clines = dset$clines[keep,]
        dset$copies = dset$copies[,keep]
        dset$eset = dset$eset[,keep]
        dset$meth = dset$meth[,keep]
        covar = 1
    }

    eset = DESeq2::estimateSizeFactors(dset$eset_raw)
    copies = dset$copies

    cnts = DESeq2::counts(eset)
    names(dimnames(cnts)) = names(dimnames(copies)) = c("gene", "cell_line")
    sfs = data.frame(cell_line=colnames(eset), sf=DESeq2::sizeFactors(eset))
    cv = data.frame(cell_line=colnames(eset), covar=covar)
    cnts = reshape2::melt(cnts, value.name="expr")
    copies2 = reshape2::melt(copies, value.name="copies")
    df = inner_join(cnts, copies2) %>%
        inner_join(sfs) %>%
        inner_join(cv) %>%
        na.omit() %>%
        mutate(eup_equiv = (copies - 2) / 2) %>%
        group_by(gene) %>%
        tidyr::nest() %>%
        tidyr::expand_grid(tibble(cna = c("amp", "del", "all")))



#    w = clustermq::workers(n_jobs = as.integer(args$cores),
#                           template = list(memory = as.integer(args$memory)))

top100 = readxl::read_xlsx("../merge_rank/rank_top/pan_rlm3.xlsx", "amp") %>%
    pull(name) %>% head(100)
genes = genes[genes %in% top100]

    ffuns = list(
        amp = function(x) { x[x < 2-et] = NA; x },
        del = function(x) { x[x > 2+et] = NA; x },
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        res = clustermq::Q(do_fit, genes=genes, n_jobs=10, memory=3072, #workers=w,
                           chunk_size=1, pkgs=c("dplyr", "rstanarm"),
                const = list(eset=eset, copies=ff(dset$copies),
                             covar=dset$clines$tcga_code, et=et)) %>%
            setNames(names(genes)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
})
