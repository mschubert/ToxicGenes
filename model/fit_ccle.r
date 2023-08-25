library(brms)
library(dplyr)
sys = import('sys')

#' Fit a negative binomial model using a stan glm
#'
#' @param dset  A data.frame with columns: expr, copies, sf, covar, eup_equiv
#' @param cna   Character string of either: "amp", "del", or "all"
#' @param mod   Precompiled brms model
#' @param et    Tolerance in ploidies to consider sample euploid
#' @param min_aneup  Minimum number of aneuploid (triploid minus et) samples
#' @return      A data.frame with fit statistics
do_fit = function(dset, cna, mod, et=0.15, min_aneup=3) {
    if (cna == "amp") {
        dset = dset %>% dplyr::filter(copies > 2-et)
    } else if (cna == "del") {
        dset = dset %>% dplyr::filter(copies < 2+et)
    }

    dset$sf = dset$sf * mean(dset$expr) # parameterize so fitted mean is constant

    n_aneup = sum(abs(dset$copies-2) > 1-et)
    if (n_aneup < min_aneup || all(dset$expr == 0))
        return(data.frame(n_aneup=n_aneup))

    init_fun = function() {
        list(b_scaling = array(runif(length(levels(dset$covar)), 0.5, 1.5)),
             b_deviation = array(0))
    }
    res = update(mod, newdata=dset, chains=4, iter=2000, init=init_fun)

    rmat = as.matrix(res)
    is_covar = grepl("covar", colnames(rmat), fixed=TRUE)
    intcp = rmat[,colnames(rmat) == "b_scaling_sf:eup_equiv" | is_covar, drop=FALSE]
    eup_dev = rmat[,grepl("eup_dev", colnames(rmat), fixed=TRUE)]
    z_comp = mean(eup_dev) / sd(eup_dev)

    tibble(estimate = mean(eup_dev) / mean(intcp),
           std.error = sd(eup_dev) / mean(intcp),
           z_comp = z_comp,
           n_aneup = n_aneup,
           eup_reads = round(mean(intcp) * mean(dset$expr)),
           n_eff = neff_ratio(res, pars="b_deviation_sf:eup_dev"),
           Rhat = rhat(res, pars="b_deviation_sf:eup_dev"),
           p.value = 2 * pnorm(abs(z_comp), lower.tail=FALSE))
}

#' Precompile BRMS model
#'
#' Using nonlinear syntax to specify coefficient bounds:
#'   https://github.com/paul-buerkner/brms/issues/1422
#'
#' @param data  A data.frame with the model data
#' @return      A brms model object
make_mod = function(data) {
    if (length(unique(data$covar)) == 1) {
        fml = bf(expr ~ 0 + scaling + deviation,
                 scaling ~ 0 + sf:eup_equiv,
                 deviation ~ 0 + sf:eup_dev,
                 nl = TRUE)
    } else {
        fml = bf(expr ~ 0 + scaling + deviation,
                 scaling ~ 0 + sf:covar:eup_equiv,
                 deviation ~ 0 + sf:eup_dev,
                 nl = TRUE)
    }

    mod = brm(fml, family=negbinomial(link="identity"),
              data = data, chains = 0, cores = 1, drop_unused_levels = FALSE,
              prior = prior(normal(0,0.5), class="b", nlpar="deviation") +
                      prior(lognormal(0,1), class="b", nlpar="scaling", lb=0))
}

#' Prepare common CCLE data.frame to serve as model input
#'
#' @param ccle_df  data.frame from ccle
#' @param tissue   character vector of tissue types, or 'pan'
prep_data = function(ccle_df, tissue) {
    if (length(tissue) == 1 && tissue == "NSCLC")
        tissue = c("LUAD", "LUSC")
    if (!identical(tissue, "pan"))
        ccle_df = ccle_df %>% filter(covar %in% tissue)

    ccle_df %>%
        mutate(eup_dev = ((copies - 2) / 2),
               eup_equiv = eup_dev + 1,
               covar = factor(covar)) %>%
        group_by(gene) %>%
            tidyr::nest() %>%
        ungroup()
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'rds', '../data/df_ccle.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('j', 'cores', 'integer', '70'),
        opt('m', 'memory', 'integer', '1024'),
        opt('o', 'outfile', 'rds', 'fit_brms/pan.rds')
    )

    cna_cmq = function(.gene, .data, cna) {
        clustermq::Q(do_fit, dset=.data,
                     const = list(cna=cna, et=cfg$euploid_tol, mod=mod),
                     pkgs = c("dplyr", "brms"),
                     n_jobs = as.integer(args$cores),
                     memory = as.integer(args$memory),
                     chunk_size = 1) %>%
            bind_rows() %>%
            mutate(gene = .gene,
                   adj.p = p.adjust(p.value, method="fdr")) %>%
            select(gene, everything()) %>%
            arrange(adj.p, p.value)
    }

    cfg = yaml::read_yaml(args$config)
    df = readRDS(args$infile) %>% prep_data(args$tissue)
    mod = make_mod(df$data[[1]])

    res = with(df, list(
        amp = cna_cmq(gene, data, "amp"),
        del = cna_cmq(gene, data, "del"),
        all = cna_cmq(gene, data, "all")
    ))

    saveRDS(res, args$outfile)
})
