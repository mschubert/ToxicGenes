library(brms)
library(dplyr)
#library(cmdstanr)
snb = import('../../ccle/fit_stan-nb')

do_fit = function(dset, et=0.15) {
    dset = dset %>% dplyr::filter(copies > 2-et)
    dset$sf = dset$sf * mean(dset$expr)
    n_aneup = sum(abs(dset$copies-2) > 1-et)

    if (length(unique(dset$covar)) == 1) {
        fml = expr ~ 0 + sf:eup_equiv + sf:eup_dev
    } else {
        fml = expr ~ 0 + sf:covar:eup_equiv + sf:eup_dev
    }
    sc = make_stancode(fml, family=negbinomial(link="identity"), data=dset,
                       prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                               prior(lognormal(0,1), class="b"))
    init_fun = function() {
        lp = grep("lprior \\+= .*b\\[[0-9]+\\]", strsplit(sc, "\\n")[[1]], value=TRUE)
        list(b = ifelse(grepl(" normal_lpdf", lp), 0, runif(length(lp), 0.5, 1.5)))
    }
    res = brm(fml, family=negbinomial(link="identity"),
              data = dset, chains=4, cores = 1, init = init_fun,
              prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                      prior(lognormal(0,1), class="b"))#,
#              backend = "cmdstanr")

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
}

ccledata = readRDS("../../data/ccle/dset.rds")
meta = ccledata$clines %>%
    select(cell_line=CCLE_ID, cline=Name) %>%
    filter(!is.na(cline)) %>% distinct()
ccle_mut = ccledata$mut %>%
    inner_join(meta) %>%
    filter(gene == "TP53",
           ! type %in% c("Silent", "3'UTR", "5'UTR", "5'Flank", "IGR", "Intron"))
ccle_df = readRDS("../../data/df_ccle.rds") %>%
    mutate(p53_mut = ifelse(cell_line %in% ccle_mut$cell_line, 1, 0)) %>%
    filter(gene == "CDKN1A")

df = snb$prep_data(ccle_df, "pan") %>% tidyr::unnest(data)
df2 = snb$prep_data(ccle_df, "BRCA") %>% tidyr::unnest(data)
df3 = snb$prep_data(ccle_df, c("LUAD", "LUSC", "SCLC")) %>% tidyr::unnest(data)

r1 = do_fit(df)
r2 = do_fit(df2)
r3 = do_fit(df3)
