library(brms)
library(dplyr)
sys = import('sys')
tcga = import('data/tcga')

#' Fit a negative binomial model using a stan glm
#'
#' @param dset  A data.frame with columns: purity, expr, copies, sf, covar, eup_equiv
#' @param cna   Character string of either: "amp", "del", or "all"
#' @param mod   A precompiled brms model (depending on fit type)
#' @param et    Tolerance in ploidies to consider sample euploid
#' @param min_aneup  Minimum number of aneuploid (triploid minus et) samples
#' @return      A data.frame with fit statistics
do_fit = function(dset, cna, mod, et=0.15, min_aneup=5) {
    stopifnot(requireNamespace("brms"))
    if (cna == "amp") {
        dset = dset %>% dplyr::filter(copies > 2-et)
    } else if (cna == "del") {
        dset = dset %>% dplyr::filter(copies < 2+et)
    }

    dset$sf = dset$sf * mean(dset$expr)

    n_aneup = sum(abs(dset$cancer_copies-2) > 1-et)
    if (n_aneup < min_aneup)
        return(data.frame(n_aneup=n_aneup))

    tryCatch({
        # stancode(mod) has eup_dev on different positions of b[i]
        init_fun = function() {
            lp = grep("lprior \\+= .*b\\[[0-9]+\\]",
                      strsplit(stancode(mod), "\\n")[[1]], value=TRUE)
            list(b = ifelse(grepl(" normal_lpdf", lp), 0, runif(length(lp), 0.5, 1.5)))
        }
        res = update(mod, newdata=dset, chains=4, iter=2000, init=init_fun)

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
               eup_reads = mean(intcp) * mean(dset$expr),
               stroma_reads = mean(stroma) * mean(dset$expr),
               n_eff = neff_ratio(res, pars=grep("eup_dev", colnames(rmat), value=TRUE)),
               Rhat = rhat(res, pars=grep("eup_dev", colnames(rmat), value=TRUE)),
               p.value = 2 * pnorm(abs(z_comp), lower.tail=FALSE))

    }, error = function(e) {
        warning(conditionMessage(e), immediate.=TRUE)
        tibble(n_aneup = n_aneup)
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
    } else {
        ex = "sf:purity:eup_dev_cancer"
        if (length(unique(data$covar)) == 1) {
            fml = expr ~ 0 + sf:purity:eup_equiv_cancer + sf:stroma + sf:purity:eup_dev_cancer
        } else if (type == "pur") {
            fml = expr ~ 0 + sf:purity:covar:eup_equiv_cancer + sf:stroma + sf:purity:eup_dev_cancer
        } else if (type == "puradj") {
            fml = expr ~ 0 + sf:purity:covar:eup_equiv_cancer + sf:covar:stroma + sf:purity:eup_dev_cancer
        }
    }

    brm(fml, family = negbinomial(link="identity"),
        data = data, chains = 0, cores = 1,
        prior = prior_string("normal(0,0.2)", coef=ex) +
                prior(lognormal(0,1), class="b"))
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
        mutate(stroma = 1 - purity,
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
        opt('t', 'tissue', 'TCGA identifier', 'COADREAD'),
        opt('y', 'type', 'naive|pur|puradj', 'pur'),
        opt('j', 'cores', 'integer', '50'),
        opt('m', 'memory', 'integer', '1024'),
        opt('o', 'outfile', 'xlsx', 'COADREAD/stan-nb_pur.xlsx')
    )

    cna_cmq = function(.gene, .data, cna) {
        clustermq::Q(do_fit, dset=.data,
                     const = list(cna=cna, mod=mod, et=cfg$euploid_tol),
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
    mod = make_mod(df$data[[1]], type=args$type)

    #TODO: remove head; for debug
    res = with(head(df), list(
        amp = cna_cmq(gene, data, "amp"),
        del = cna_cmq(gene, data, "del"),
        all = cna_cmq(gene, data, "all")
    ))

    writexl::write_xlsx(res, args$outfile)
})
