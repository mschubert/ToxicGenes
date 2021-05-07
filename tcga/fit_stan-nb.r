library(brms)
library(dplyr)
sys = import('sys')
tcga = import('data/tcga')

#' Fit a negative binomial model using a stan glm
#'
#' @param df   A data.frame with columns: purity, expr, copies, sf, covar, eup_equiv
#' @param cna  Character string of either: "amp", "del", or "all"
#' @param mod  A precompiled brms model (depending on fit type)
#' @param et   Tolerance in ploidies to consider sample euploid
#' @param min_aneup  Minimum number of aneuploid (triploid minus et) samples
#' @return     A data.frame with fit statistics
do_fit = function(df, cna, mod, et=0.15, min_aneup=5) {
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
        init_vars = with(prior_summary(mod), coef[class == "b" & coef != ""])
        init_fun = function() list(b=runif(length(init_vars), 0.5, 1.5))

        res = update(mod, newdata=df, chains=4, iter=2000, init=init_fun)

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
        tibble(estimate = NA)
    })
}

#' Precompile BRMS model
#'
#' @param data  A data.frame with the model data
#' @param type  Character string of model type, i.e. 'naive' or 'pur'
#' @return      A brms model object
make_mod = function(data, type="naive") {
    if (type == "naive") {
        ex = "sf:eup_dev"
        if (length(unique(data$covar)) == 1) {
            fml = expr ~ 0 + sf:eup_equiv + sf:eup_dev
        } else {
            fml = expr ~ 0 + sf:covar:eup_equiv + sf:eup_dev
        }
    } else if (type == "pur") {
        ex = "sf:purity:eup_dev_cancer"
        if (length(unique(data$covar)) == 1) {
            fml = expr ~ 0 + sf:purity:eup_equiv_cancer + sf:stroma + sf:purity:eup_dev_cancer
        } else {
            fml = expr ~ 0 + sf:purity:covar:eup_equiv_cancer + sf:stroma + sf:purity:eup_dev_cancer
        }
    } else
        stop("invalid model 'type'") #todo: stroma expr per type?

    brm(fml, family = negbinomial(link="identity"),
        data = data, chains = 0, cores = 1,
        prior = prior_string("normal(0,0.2)", coef=ex) +
                prior(gamma(1.5,1), class="b"))
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

    cna_cmq = function(data, cna, mod) {
        clustermq::Q(do_fit, df=data,
                     const = list(cna=cna, mod=mod, et=cfg$euploid_tol),
                     pkgs = c("dplyr", "brms"),
                     n_jobs = as.integer(args$cores),
                     memory = as.integer(args$memory),
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
        ungroup() %>%
        mutate(amp = cna_cmq(data, "amp", make_mod(data[[1]], args$type)))
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
