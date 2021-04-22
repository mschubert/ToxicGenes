library(dplyr)
sys = import('sys')

#' Fit a negative binomial model using a stan glm
#'
#' @param df   A data.frame with columns: expr, copies, sf, covar, eup_equiv
#' @param cna  Character string of either: "amp", "del", or "all"
#' @param et   Tolerance in ploidies to consider sample euploid
#' @param timeout  Number of seconds that a fit can take before returnin NA
#' @return     A data.frame with fit statistics
do_fit = function(df, cna, et=0.15, timeout=1800) {
    stopifnot(requireNamespace("rstanarm"))

    if (cna == "amp") {
        df = df %>% filter(copies > 2-et)
    } else if (cna == "del") {
        df = df %>% filter(copies < 2+et)
    }

    if (length(unique(df$covar)) == 1) {
        full = expr ~ sf + eup_dev:sf
    } else {
        full = expr ~ sf + covar:sf + eup_dev:sf
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
        intcp = rmat[,colnames(rmat) == "sf" | is_covar, drop=FALSE] + rmat[,"(Intercept)"]
        sd_intcp = mean(apply(intcp, 2, sd))
        eup_eq = rmat[,grepl("eup_dev", colnames(rmat), fixed=TRUE)] # term can be: eup_dev:sf, sf:eup_dev
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
        data.frame(estimate=NA)
    })
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'rds', '../data/df_ccle.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('j', 'cores', 'integer', '20'),
        opt('m', 'memory', 'integer', '1024'),
        opt('o', 'outfile', 'xlsx', 'pan/stan-nb.xlsx')
    )

    cna_cmq = function(data, cna) {
        clustermq::Q(do_fit, df=data,
                     const = list(cna=cna, et=cfg$euploid_tol),
                     pkgs = c("dplyr", "rstanarm"),
                     n_jobs = as.integer(args$cores),
                     memory = as.integer(args$memory),
#                     max_calls_worker = 142, # 72h job / 0.5h max run - 1h
                     chunk_size = 1)
    }

    cfg = yaml::read_yaml(args$config)
    df = readRDS(args$infile)
    if (args$tissue == "NSCLC") {
        args$tissue = c("LUAD", "LUSC")
    }
    if (!identical(args$tissue, "pan")) {
        df = df %>% filter(covar %in% args$tissue)
    }
    df = df %>%
        mutate(eup_dev = eup_equiv - 1) %>%
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
