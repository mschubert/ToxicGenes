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
    stopifnot(requireNamespace("brms"))
    if (cna == "amp") {
        dset = dset %>% dplyr::filter(copies > 2-et)
    } else if (cna == "del") {
        dset = dset %>% dplyr::filter(copies < 2+et)
    }

    dset$sf = dset$sf * mean(dset$expr) # parameterize so fitted mean is constant

    n_aneup = sum(abs(dset$copies-2) > 1-et)
    if (n_aneup < min_aneup)
        return(data.frame(n_aneup=n_aneup))

    tryCatch({
        # stancode(mod) has eup_dev first and then covar intercepts
        init_fun = function() list(b=c(0, runif(length(unique(dset$covar)), 0.5, 1.5)))
        res = update(mod, newdata=dset, chains=4, iter=2000, init=init_fun)

        rmat = as.matrix(res)
        is_covar = grepl("covar", colnames(rmat), fixed=TRUE)
        intcp = rmat[,colnames(rmat) == "b_sf:eup_equiv" | is_covar, drop=FALSE]
        eup_dev = rmat[,grepl("eup_dev", colnames(rmat), fixed=TRUE)]
        z_comp = mean(eup_dev) / sd(eup_dev)

        tibble(estimate = mean(eup_dev) / mean(intcp),
               z_comp = z_comp,
               n_aneup = n_aneup,
               n_genes = 1,
               eup_reads = mean(intcp) * mean(dset$expr),
               n_eff = neff_ratio(res, pars="b_sf:eup_dev"),
               Rhat = rhat(res, pars="b_sf:eup_dev"),
               p.value = 2 * pnorm(abs(z_comp), lower.tail=FALSE))

    }, error = function(e) {
        warning(conditionMessage(e), immediate.=TRUE)
        tibble(n_aneup=n_aneup)
    })
}

#' Precompile BRMS model
#'
#' @param data  A data.frame with the model data
#' @return      A brms model object
make_mod = function(data) {
    if (length(unique(df$data[[1]]$covar)) == 1) {
        fml = expr ~ 0 + sf:eup_equiv + sf:eup_dev
    } else {
        fml = expr ~ 0 + sf:covar:eup_equiv + sf:eup_dev
    }

    mod = brm(fml, family=negbinomial(link="identity"),
              data = data, chains = 0, cores = 1,
              prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                      prior(lognormal(0,1), class="b"))
}

#' Prepare common CCLE data.frame to serve as model input
#'
#' @param ccle_df  data.frame from ccle
#' @param tissue   character vector of tissue types, or 'pan'
prep_data = function(ccle_df, tissue) {
    if (tissue == "NSCLC")
        tissue = c("LUAD", "LUSC")
    if (!identical(args$tissue, "pan"))
        df = df %>% filter(covar %in% tissue)

    ccle_df %>%
        mutate(eup_dev = ((copies - 2) / 2),
               eup_equiv = eup_dev + 1) %>%
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
        opt('o', 'outfile', 'xlsx', 'pan/stan-nb.xlsx')
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

    writexl::write_xlsx(res, args$outfile)
})
