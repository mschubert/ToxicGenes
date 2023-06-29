library(brms)
library(dplyr)
snb = import('../../ccle/fit_stan-nb')

df = readRDS("../../data/df_ccle.rds") %>%
    snb$prep_data("pan") %>%
    filter(gene == "CDKN1A") %>%
    tidyr::unnest(data)

mod = snb$make_mod(df)

res1 = snb$do_fit(df, "amp", mod)
