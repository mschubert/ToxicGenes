library(brms)
library(dplyr)
sys = import('sys')

#' Fit a negative binomial model using a stan glm
#'
#' @param df   A data.frame with columns: expr, copies, sf, covar, eup_equiv
#' @param cna  Character string of either: "amp", "del", or "all"
#' @param mod  Precompiled brms model
#' @param et   Tolerance in ploidies to consider sample euploid
#' @param timeout  Number of seconds that a fit can take before returnin NA
#' @return     A data.frame with fit statistics
do_fit = function(df, cna, mod, mod_covar, et=0.15, min_aneup=3, timeout=1800) {
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

    tryCatch({
        setTimeLimit(elapsed=timeout, transient=TRUE)

        init_vars = sum(grepl("eup_equiv", prior_summary(mod)$coef)) + 1
        init_fun = function() list(b=runif(init_vars, 0.5, 2))

        res = update(mod, newdata=df, chains=4, iter=2000, seed=2380719, init=init_fun)

        rmat = as.matrix(res)
        is_covar = grepl("covar", colnames(rmat), fixed=TRUE)
        intcp = rmat[,colnames(rmat) == "b_sf:eup_equiv" | is_covar, drop=FALSE]
        eup_eq = rmat[,grepl("eup_dev", colnames(rmat), fixed=TRUE)]
        z_comp = mean(eup_eq) / sd(eup_eq)

        tibble(estimate = mean(eup_eq) / mean(intcp),
               z_comp = z_comp,
               n_aneup = n_aneup,
               n_genes = 1,
               eup_reads = mean(intcp) * mean(df$expr),
               slope_diff = (mean(eup_eq) - mean(intcp))/2 * mean(df$expr),
               p.value = 2 * pnorm(abs(z_comp), lower.tail=FALSE))

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
                     const = list(cna=cna, et=cfg$euploid_tol, mod=mod),
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
    if (length(unique(df$data[[1]]$covar)) == 1) {
        mod = brm(expr ~ 0 + sf:eup_equiv + sf:eup_dev, family=negbinomial(link="identity"),
                  data=df$data[[1]], chains=0, cores=1,
                  prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                          prior(normal(1,5), class="b"))
    } else {
        mod = brm(expr ~ 0 + sf:covar:eup_equiv + sf:eup_dev, family=negbinomial(link="identity"),
                  data=df$data[[1]], chains=0, cores=1,
                  prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                          prior(normal(1,5), class="b"))
    }
    df = df %>%
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
