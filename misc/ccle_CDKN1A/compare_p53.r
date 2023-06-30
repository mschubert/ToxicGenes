library(brms)
library(dplyr)
library(cmdstanr)
snb = import('../../ccle/fit_stan-nb')

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

fml = expr ~ 0 + sf:covar:eup_equiv + sf:eup_dev

if2 = function() {
    list(b=c(runif(1, -0.1, 0.1), runif(length(unique(df$covar)), 0.5, 1.5)))
}

mod = brm(fml, family=negbinomial(link="identity"),
          data = df, chains=0, cores = 1, init = if2,
          prior = prior(normal(0,0.5), coef="sf:eup_dev") +
                  prior(lognormal(0,1), class="b"),
          backend = "cmdstanr")
