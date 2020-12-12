library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

#' Fit a negative binomial model using a stan glm
#'
#' @param df    A data.frame with columns: purity, expr, copies, sf, covar, eup_equiv
#' @param cna   Character string of either: "amp", "del", or "all"
#' @param type  Character string of either: "naive", "pur"
#' @param et    Tolerance in ploidies to consider sample euploid
#' @param timeout  Number of seconds that a fit can take before returnin NA
#' @return     A data.frame with fit statistics
do_fit = function(df, cna, type="pur", et=0.15, timeout=3600) {
    stopifnot(requireNamespace("rstanarm"))
    stopifnot(requireNamespace("statip"))

    if (cna == "amp") {
        df = df %>% filter(copies > 2-et)
    } else if (cna == "del") {
        df = df %>% filter(copies < 2+et)
    }

    if (length(unique(df$covar)) == 1) {
        full = switch(type,
            naive = expr ~ eup_equiv,
            pur = expr ~ purity + eup_equiv
#            puradj = expr ~ purity + eup_equiv
        )
    } else {
        full = switch(type,
            naive = expr ~ covar + eup_equiv,
            pur = expr ~ covar + purity + eup_equiv
#            puradj = expr ~ covar + purity + eup_equiv
        )
    }

    tryCatch({
        setTimeLimit(elapsed=timeout, transient=TRUE)
        res = rstanarm::stan_glm(full, data=df, offset=sizeFactor,
                                 family = neg_binomial_2(link="identity"),
                                 prior = normal(mean(df$expr), sd(df$expr)),
                                 prior_intercept = normal(mean(df$expr), sd(df$expr)),
                                 seed = 2380719)

        rmat = as.matrix(res)
        hdist = statip::hellinger(rmat[,"(Intercept)"], rmat[,"eup_equiv"])

        pseudo_p = pnorm(abs((mean(rmat[,"(Intercept)"]) - mean(rmat[,"eup_equiv"])) /
                             sd(rmat[,"(Intercept)"])), lower.tail=F)

        tibble(estimate = mean(rmat[,"eup_equiv"]) / mean(rmat[,"(Intercept)"]) - 1, # pct_comp
               n_aneup = sum(abs(df$copies-2) > 1-et),
               n_genes = 1,
               eup_reads = mean(rmat[,"(Intercept)"]),
               slope_diff = mean(rmat[,"(Intercept)"]) - mean(rmat[,"eup_equiv"]),
               rsq = hdist, # not rsq, but [0,1]
               p.value = pseudo_p)

    }, error = function(e) {
        warning(conditionMessage(e), immediate.=TRUE)
        data.frame(estimate = NA)
    })
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'rds', '../data/df_tcga.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'LUAD'),
        opt('y', 'type', 'naive|pur|puradj', 'naive'),
        opt('j', 'cores', 'integer', '1000'),
        opt('m', 'memory', 'integer', '1024'),
        opt('o', 'outfile', 'xlsx', 'LUAD/genes.xlsx')
    )

    cna_cmq = function(data, cna) {
        clustermq::Q(do_fit, df=data,
                     const = list(cna=cna, type=args$type, et=cfg$euploid_tol),
                     pkgs = c("dplyr", "rstanarm"),
                     n_jobs = as.integer(args$cores),
                     memory = as.integer(args$memory),
                     chunk_size = 1)
    }

    cfg = yaml::read_yaml(args$config)
    if (args$tissue == "pan")
        args$tissue = tcga$cohorts()

    df = readRDS(args$infile) %>%
        filter(covar %in% args$tissue) %>%
        group_by(gene) %>%
        tidyr::nest() %>%
        ungroup() %>%
        mutate(amp = cna_cmq(data, "amp"))
#               del = cna_cmq(data, "del"),
#               all = cna_cmq(data, "all"))

    saveRDS(df, file=args$outfile)

#    res = df %>%
#        select(-data) %>%
#        tidyr::unnest("res") %>%
#        split(.$cna) %>%
#        lapply(. %>% mutate(adj.p = p.adjust(p.value, method="fdr")) %>% arrange(adj.p, p.value))

#    writexl::write_xlsx(res, args$outfile)
})
