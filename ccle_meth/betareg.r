library(dplyr)
library(betareg)
sys = import('sys')

fit_beta = function(df) {
    df = na.omit(df[c("cohort", "meth", "copies")])
    df$meth = pmin(df$meth, 1 - 1e-5)
    df$meth = pmax(df$meth, 1e-5)

    tryCatch({
            betareg(meth ~ cohort + copies, data=df) %>%
                broom::tidy() %>%
                filter(term == "copies") %>%
                select(-component, -term)
        }, error = function(e) data.frame(estimate=NA))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
    opt('c', 'cohort', 'pan|TCGA cohort', 'pan'),
    opt('e', 'euploid_tol', 'numeric copy dev from euploid', '0.15'),
    opt('o', 'outfile', 'rds', 'betareg.rds')
)

#TODO:
#  volcano plot
# gene sets + plot (different script?)
#  overlap ccle stan-nb compensation

et = as.numeric(args$euploid_tol)
ccle = readRDS(args$infile)

names(dimnames(ccle$copies)) = names(dimnames(ccle$meth)) = c("gene", "CCLE_ID")
cdf = ccle$clines %>%
    select(CCLE_ID, Name, Site_Primary, cohort=tcga_code) %>%
    left_join(reshape2::melt(ccle$copies, value.name="copies")) %>%
    left_join(reshape2::melt(ccle$meth, value.name="meth")) %>%
    filter(copies > 2-et) %>%
    group_by(gene) %>%
    tidyr::nest()

res = cdf %>%
    rowwise() %>%
        mutate(res = list(fit_beta(data))) %>%
    ungroup() %>%
    tidyr::unnest("res") %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)
#res2 = res %>% filter(gene != "LDLR")

saveRDS(res, file=args$outfile)
