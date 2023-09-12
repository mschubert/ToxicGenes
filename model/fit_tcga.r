library(brms)
library(dplyr)
sys = import('sys')
tcga = import('data/tcga')

#' Fit a negative binomial model using a stan glm
#'
#' @param dset  A data.frame with columns: purity, expr, copies, sf, covar, eup_equiv
#' @param cna   Character string of either: "amp", "del", or "all"
#' @param mod   A precompiled brms model (depending on fit type)
#' @param type  The fit type, ie. 'naive|pur|puradj'
#' @param et    Tolerance in ploidies to consider sample euploid
#' @param min_aneup  Minimum number of aneuploid (triploid minus et) samples
#' @return      A data.frame with fit statistics
do_fit = function(dset, cna, mod, type, et=0.15, min_aneup=5) {
    if (cna == "amp") {
        dset = dset %>% dplyr::filter(copies > 2-et)
    } else if (cna == "del") {
        dset = dset %>% dplyr::filter(copies < 2+et)
    }

    dset$sf = dset$sf * mean(dset$expr)

    n_aneup = sum(abs(dset$cancer_copies-2) > 1-et)
    if (n_aneup < min_aneup || all(dset$expr == 0))
        return(data.frame(n_aneup=n_aneup))

    init_fun = function() {
        n_noncancer = switch(type, naive=0, pur=1, puradj=length(levels(dset$covar)))
        list(b_scaling = array(runif(length(levels(dset$covar)), 0.5, 1.5)),
             b_noncancer = array(runif(n_noncancer, 0.5, 1.5)),
             b_deviation = array(0))
    }
    res = update(mod, newdata=dset, chains=4, iter=2000, init=init_fun)

    rmat = as.matrix(res)
    is_covar = grepl("covar", colnames(rmat))
    intcp = rmat[,grepl("eup_equiv", colnames(rmat)) | is_covar, drop=FALSE]
    eup_dev = rmat[,grepl("eup_dev", colnames(rmat), fixed=TRUE)]
    stroma = rmat[,grepl("stroma", colnames(rmat), fixed=TRUE)]
    z_comp = mean(eup_dev) / sd(eup_dev)

    tibble(estimate = mean(eup_dev) / mean(intcp),
           std.error = sd(eup_dev) / mean(intcp),
           z_comp = z_comp,
           n_aneup = n_aneup,
           eup_reads = round(mean(intcp) * mean(dset$expr)),
           stroma_reads = round(mean(stroma) * mean(dset$expr)),
           n_eff = neff_ratio(res, pars=grep("eup_dev", colnames(rmat), value=TRUE)),
           Rhat = rhat(res, pars=grep("eup_dev", colnames(rmat), value=TRUE)),
           p.value = 2 * pnorm(abs(z_comp), lower.tail=FALSE))
}

#' Precompile BRMS model
#'
#' @param data  A data.frame with the model data
#' @param type  Character string of model type, i.e. 'naive' or 'pur'
#' @return      A brms model object
make_mod = function(data, type="naive") {
    if (type == "naive") {
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
        add_stroma_prior = c()
    } else {
        if (length(unique(data$covar)) == 1) {
            fml = bf(expr ~ 0 + scaling + noncancer + deviation,
                     scaling ~ 0 + sf:purity:eup_equiv_cancer,
                     noncancer ~ 0 + sf:stroma,
                     deviation ~ 0 + sf:purity:eup_dev_cancer,
                     nl = TRUE)
        } else if (type == "pur") {
            fml = bf(expr ~ 0 + scaling + noncancer + deviation,
                     scaling ~ 0 + sf:purity:covar:eup_equiv_cancer,
                     noncancer ~ 0 + sf:stroma,
                     deviation ~ 0 + sf:purity:eup_dev_cancer,
                     nl = TRUE)
        } else if (type == "puradj") {
            fml = bf(expr ~ 0 + scaling + noncancer + deviation,
                     scaling ~ 0 + sf:purity:covar:eup_equiv_cancer,
                     noncancer ~ 0 + sf:covar:stroma,
                     deviation ~ 0 + sf:purity:eup_dev_cancer,
                     nl = TRUE)
        }
        add_stroma_prior = prior(lognormal(0,1), class="b", nlpar="noncancer", lb=0)
    }

    brm(fml, family = negbinomial(link="identity"),
        data = data, chains = 0, cores = 1, drop_unused_levels = FALSE,
        prior = prior(normal(0,0.2), class="b", nlpar="deviation") +
                prior(lognormal(0,1), class="b", nlpar="scaling", lb=0) +
                add_stroma_prior)
}

#' Prepare common TCGA data.frame to serve as model input
#'
#' @param tcga_df
#' @param tissue
prep_data = function(tcga_df, tissue) {
    tissue = switch(tissue,
        pan = tcga$cohorts(),
        COADREAD = c("COAD", "READ"),
        NSCLC = c("LUAD", "LUSC"),
        tissue
    )

    tcga_df %>%
        filter(covar %in% tissue) %>%
        mutate(covar = factor(covar),
               stroma = 1 - purity,
               eup_dev = (copies - 2) / 2,
               eup_equiv = eup_dev + 1,
               eup_dev_cancer = (cancer_copies - 2) / 2,
               eup_equiv_cancer = eup_dev_cancer + 1) %>%
        group_by(gene) %>%
            tidyr::nest() %>%
        ungroup()
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'rds', '../data/df_tcga.rds'),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('p', 'cna', 'amp|del|all', 'amp'),
        opt('y', 'type', 'naive|pur|puradj', 'puradj'),
        opt('j', 'cores', 'integer', '70'),
        opt('m', 'memory', 'integer', '1536'),
        opt('o', 'outfile', 'xlsx', 'fit_tcga_puradj-amp/pan.rds')
    )

    cfg = yaml::read_yaml(args$config)
    df = readRDS(args$infile) %>% prep_data(args$tissue)
    mod = make_mod(df$data[[1]], type=args$type)
    const = c(list(mod=mod, et=cfg$euploid_tol), args[c("cna", "type")])

    res = clustermq::Q(do_fit, dset = df$data,
                       const = const, pkgs = c("dplyr", "brms"),
                       n_jobs = as.integer(args$cores),
                       memory = as.integer(args$memory),
                       chunk_size = 1) %>%
        setNames(as.character(df$gene)) %>%
        bind_rows(.id="gene") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)

    saveRDS(res, args$outfile)
})
