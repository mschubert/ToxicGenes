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
do_fit = function(df, cna, type="pur", et=0.15, timeout=7200) {
    stopifnot(requireNamespace("rstanarm"))

    if (cna == "amp") {
        df = df %>% filter(copies > 2-et)
    } else if (cna == "del") {
        df = df %>% filter(copies < 2+et)
    }
    df$sf = df$sizeFactor
    df$stroma = 1 - df$purity

    if (length(unique(df$covar)) == 1) {
        full = switch(type,
            naive = expr ~ eup_equiv:sf,
            pur = expr ~ purity:sf + eup_equiv:sf,
            puradj = expr ~ stroma:sf + purity:eup_equiv:sf
        )
    } else {
        full = switch(type,
            naive = expr ~ covar:sf + eup_equiv:sf,
            pur = expr ~ covar:sf + purity:sf + eup_equiv:sf,
            puradj = expr ~ covar:sf + stroma:sf + purity:eup_equiv:sf
#            puradj = expr ~ covar:sf + covar:purity:sf + purity:eup_equiv:sf
        )
    }

    tryCatch({
        setTimeLimit(elapsed=timeout, transient=TRUE)
        res = rstanarm::stan_glm(full, data=df,
                                 family = neg_binomial_2(link="identity"),
                                 prior = normal(mean(df$expr), sd(df$expr)),
                                 prior_intercept = normal(mean(df$expr), sd(df$expr)),
                                 seed = 2380719)

        rmat = as.matrix(res)
        is_covar = grepl("covar", colnames(rmat))
        rmat[,is_covar] = rmat[,is_covar] + rmat[,"(Intercept)"]
        intcp = rmat[colnames(rmat) == "(Intercept)" | is_covar,, drop=FALSE]
        sd_intcp = mean(apply(intcp, 2, sd))
        eup_eq = rmat[,grepl("eup_equiv", colnames(rmat), fixed=TRUE)]
        pseudo_p = pnorm(abs((mean(intcp) - mean(eup_eq)) / sd_intcp), lower.tail=F)

        tibble(estimate = mean(eup_eq) / mean(intcp) - 1, # pct_comp
               n_aneup = sum(abs(df$copies-2) > 1-et),
               n_genes = 1,
               eup_reads = mean(intcp),
               slope_diff = mean(eup_eq) - mean(intcp),
               cv_intcp = sd_intcp / mean(intcp),
               cv_copy = sd(eup_eq) / mean(eup_eq),
#               rsq = hdist, # not rsq, but [0,1]
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
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('y', 'type', 'naive|pur|puradj', 'pur'),
        opt('j', 'cores', 'integer', '1000'),
        opt('m', 'memory', 'integer', '1024'),
        opt('o', 'outfile', 'xlsx', 'pan/stan-nb_pur.xlsx')
    )

    cna_cmq = function(data, cna) {
        clustermq::Q(do_fit, df=data,
                     const = list(cna=cna, type=args$type, et=cfg$euploid_tol),
                     pkgs = c("dplyr", "rstanarm"),
                     n_jobs = as.integer(args$cores),
                     memory = as.integer(args$memory),
                     max_calls_worker = 35, # 72h job / 2h max run - 1
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

    res = df %>%
        select(-data) %>%
        tidyr::unnest("amp") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
    saveRDS(res, file=sub("\\.xlsx$", ".rds", args$outfile))

#    res = df %>%
#        select(-data) %>%
#        tidyr::unnest("res") %>%
#        split(.$cna) %>%
#        lapply(. %>% mutate(adj.p = p.adjust(p.value, method="fdr")) %>% arrange(adj.p, p.value))

    writexl::write_xlsx(res, args$outfile)
})
