library(brms)
library(dplyr)
sys = import('sys')

#' Fit a negative binomial model using a stan glm
#'
#' @param df    A data.frame with columns: expr, copies, sf, covar, eup_equiv
#' @param cna   Character string of either: "amp", "del", or "all"
#' @param mods  List of precompiled brms models
#' @param et    Tolerance in ploidies to consider sample euploid
#' @param timeout  Number of seconds that a fit can take before returnin NA
#' @return      A data.frame with fit statistics
do_fit = function(df, cna, mods, mod_covar, et=0.15, min_aneup=3, timeout=1800) {
    stopifnot(requireNamespace("brms"))

    if (cna == "amp") {
        df = df %>% filter(copies > 2-et)
    } else if (cna == "del") {
        df = df %>% filter(copies < 2+et)
    }

    df$sf = df$sf * mean(df$expr) # parameterize so prior SD is constant

    n_aneup = sum(abs(df$copies-2) > 1-et)
    if (n_aneup < min_aneup)
        return(data.frame(n_aneup=n_aneup))

    if (length(unique(df$covar)) == 1) {
        mod = mods$simple
    } else {
        mod = mods$covar
    }

    tryCatch({
        setTimeLimit(elapsed=timeout, transient=TRUE)
        res = update(mod, newdata=df, chains=4, iter=2000, seed=2380719)

        rmat = as.matrix(res)
        is_covar = grepl("covar", colnames(rmat))
#        intcp = rmat[,colnames(rmat) == "b_sf" | is_covar, drop=FALSE]
        intcp = rmat[,colnames(rmat) == "b_eup_equiv:sf" | is_covar, drop=FALSE]
        sd_intcp = mean(apply(intcp, 2, sd))
        eup_eq = rmat[,grepl("eup_dev", colnames(rmat), fixed=TRUE)] # term can be: eup_dev:sf, sf:eup_dev
        pseudo_p = pnorm(abs((mean(intcp) - mean(eup_eq)) / sd_intcp), lower.tail=F)

        tibble(estimate = mean(eup_eq) / mean(intcp) - 1, # pct_comp FIXME: why -1?
               n_aneup = n_aneup,
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
                     const = list(cna=cna, et=cfg$euploid_tol, mods=mods),
                     pkgs = c("dplyr", "brms"),
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
        group_by(gene) %>%
            tidyr::nest() %>%
        ungroup()
    testd = df$data[[1]] %>% mutate(covar = sample(letters[1:2], nrow(df$data[[1]]), replace=TRUE))
    mods = list(
        simple = brm(expr ~ 0 + eup_equiv:sf + eup_dev:sf, family=negbinomial(link="identity"),
                     data=testd, chains=1, iter=1, control=list(adapt_delta=0.99),
                     prior = prior(normal(0,0.15), coef="sf:eup_dev")),
        covar = brm(expr ~ 0 + sf + covar:sf + eup_dev:sf, family=negbinomial(link="identity"), #FIXME:
                    data=testd, chains=1, iter=1, control=list(adapt_delta=0.99))
    )
    df = df[sample(seq_len(nrow(df)), 5000),] %>%
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
