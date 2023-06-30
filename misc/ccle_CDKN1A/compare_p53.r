library(brms)
library(dplyr)
library(cmdstanr)
snb = import('../../ccle/fit_stan-nb')

df = readRDS("../../data/df_ccle.rds") %>%
    snb$prep_data("pan") %>%
    filter(gene == "CDKN1A") %>%
    tidyr::unnest(data)

#mod = make_mod(df)

#res1 = do_fit(df, "amp", mod)

fml = expr ~ 0 + sf:covar:eup_equiv + sf:eup_dev

#init_fun = function() {
#    lp = grep("lprior \\+= .*b\\[[0-9]+\\]",
#              strsplit(stancode(mod), "\\n")[[1]], value=TRUE)
#    list(b = ifelse(grepl(" normal_lpdf", lp), 0, runif(length(lp), 0.5, 1.5)))
#}

if2 = function() {
    list(b=c(runif(1, -0.1, 0.1), runif(length(unique(df$covar)), 0.5, 1.5)))
}

mod2 = brm(fml, family=negbinomial(link="identity"),
          data = df, chains=4, cores = 1, init = if2,
          prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                  prior(lognormal(0,1), class="b"),
          backend = "cmdstanr")
