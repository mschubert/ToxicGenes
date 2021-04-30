library(brms)
library(dplyr)
sys = import('sys')
tcga = import('data/tcga')

#' Fit a negative binomial model using a stan glm
#'
#' @param df    A data.frame with columns: purity, expr, copies, sf, covar, eup_equiv
#' @param cna   Character string of either: "amp", "del", or "all"
#' @param mods  A named list of precompiled brms models (names: "naive", "pur")
#' @param type  Character string of either: "naive", "pur"
#' @param et    Tolerance in ploidies to consider sample euploid
#' @param timeout  Number of seconds that a fit can take before returnin NA
#' @return     A data.frame with fit statistics
do_fit = function(df, cna, mods, type="pur", et=0.15, min_aneup=5, timeout=7200) {
    stopifnot(requireNamespace("brms"))
    if (cna == "amp") {
        df = df %>% dplyr::filter(copies > 2-et)
    } else if (cna == "del") {
        df = df %>% dplyr::filter(copies < 2+et)
    }

    df$sf = df$sf * mean(df$expr)

    n_aneup = sum(abs(df$cancer_copies-2) > 1-et)
    if (n_aneup < min_aneup)
        return(data.frame(n_aneup=n_aneup))

    tryCatch({
        setTimeLimit(elapsed=timeout, transient=TRUE)

        init_vars = with(prior_summary(mods[[type]]), coef[class == "b" & coef != ""])
        init_fun = function() list(b=runif(length(init_vars), 0.5, 1.5))

        res = update(mods[[type]], newdata=df, chains=4, iter=2000, init=init_fun)

        rmat = as.matrix(res)
        is_covar = grepl("covar", colnames(rmat))
        intcp = rmat[,grepl("eup_equiv", colnames(rmat)) | is_covar, drop=FALSE]
        eup_dev = rmat[,grepl("eup_dev", colnames(rmat), fixed=TRUE)]
        stroma = rmat[,grepl("stroma", colnames(rmat), fixed=TRUE)]
        z_comp = mean(eup_dev) / sd(eup_dev)

        tibble(estimate = mean(eup_dev) / mean(intcp),
               z_comp = z_comp,
               n_aneup = n_aneup,
               n_genes = 1,
               eup_reads = mean(intcp) * mean(df$expr),
               stroma_reads = mean(stroma) * mean(df$expr),
               p.value = 2 * pnorm(abs(z_comp), lower.tail=FALSE))

    }, error = function(e) {
        warning(conditionMessage(e), immediate.=TRUE)
        data.frame(estimate = NA)
    })
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'rds', '../data/df_tcga.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'COADREAD'),
        opt('y', 'type', 'naive|pur|puradj', 'pur'),
        opt('j', 'cores', 'integer', '50'),
        opt('m', 'memory', 'integer', '1024'),
        opt('o', 'outfile', 'xlsx', 'COADREAD/stan-nb_pur.xlsx')
    )

    cna_cmq = function(data, cna) {
        to = 3 # hours
        clustermq::Q(do_fit, df=data,
                     const = list(cna=cna, mods=mods, type=args$type,
                                  et=cfg$euploid_tol, timeout=round(to*3600)),
                     pkgs = c("dplyr", "brms"),
                     n_jobs = as.integer(args$cores),
                     memory = as.integer(args$memory),
#                     max_calls_worker = 25, #round(72 / to) - 1, # 72h job / 2h max run - 1
                     chunk_size = 1)
    }

    cfg = yaml::read_yaml(args$config)
    if (args$tissue == "pan") {
        args$tissue = tcga$cohorts()
    } else if (args$tissue == "COADREAD") {
        args$tissue = c("COAD", "READ")
    } else if (args$tissue == "NSCLC") {
        args$tissue = c("LUAD", "LUSC")
    }

    df = readRDS(args$infile) %>%
        filter(covar %in% args$tissue) %>%
        mutate(stroma = 1 - purity) %>%
        group_by(gene) %>%
            tidyr::nest() %>%
        ungroup()
    if (length(unique(df$data[[1]]$covar)) == 1) {
        mods = list(
            naive = brm(expr ~ 0 + sf:eup_equiv + sf:eup_dev,
                        family = negbinomial(link="identity"),
                        data = df$data[[1]], chains = 0, cores = 1,
                        prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                                prior(gamma(1.5,1), class="b")),
            pur = brm(expr ~ 0 + sf:purity:eup_equiv_cancer + sf:stroma + sf:eup_dev_cancer,
                      family = negbinomial(link="identity"),
                      data = df$data[[1]], chains = 0, cores = 1,
                      prior = prior(normal(0,0.5), coef="sf:eup_dev_cancer") +
                              prior(gamma(1.5,1), class="b"))
        )
    } else {
        mods = list(
            naive = brm(expr ~ 0 + sf:covar:eup_equiv + sf:eup_dev,
                        family = negbinomial(link="identity"),
                        data = df$data[[1]], chains = 0, cores = 1,
                        prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                                prior(gamma(1.5,1), class="b")),
            pur = brm(expr ~ 0 + sf:purity:covar:eup_equiv_cancer + sf:stroma + sf:eup_dev_cancer,
                      family = negbinomial(link="identity"),
                      data = df$data[[1]], chains = 0, cores = 1,
                      prior = prior(normal(0,0.5), coef="sf:eup_dev_cancer") +
                              prior(gamma(1.5,1), class="b")) #todo: stroma expr per type?
        )
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
