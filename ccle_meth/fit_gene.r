library(dplyr)
library(betareg)
sys = import('sys')
plt = import('plot')

fit_beta = function(df) {
    df = na.omit(df[c("cohort", "meth", "avg_meth_sample", "copies")])
    df$meth = pmin(df$meth, 1 - 1e-5)
    df$meth = pmax(df$meth, 1e-5)

    if (length(unique(df$cohort)) == 1)
        fml = meth ~ avg_meth_sample + copies
    else
        fml = meth ~ cohort + avg_meth_sample + copies

    tryCatch({
            betareg(fml, data=df) %>%
                broom::tidy() %>%
                filter(term == "copies") %>%
                select(-component, -term)
        }, error = function(e) data.frame(estimate=NA))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
    opt('t', 'tissue', 'pan|TCGA cohort', 'pan'),
    opt('e', 'euploid_tol', 'numeric copy dev from euploid', '0.15'),
    opt('c', 'cores', 'parallel cores', '0'),
    opt('m', 'memory', 'mem per core', '1024'),
    opt('o', 'outfile', 'rds', 'pan/betareg.rds'),
    opt('p', 'plotfile', 'pdf', 'pan/betareg.pdf')
)

et = as.numeric(args$euploid_tol)
ccle = readRDS(args$infile)
if (args$tissue == "pan")
    args$tissue = unique(na.omit(ccle$clines$tcga_code))

names(dimnames(ccle$copies)) = names(dimnames(ccle$meth)) = c("gene", "CCLE_ID")
avg_meth_sample = tibble(CCLE_ID = colnames(ccle$meth),
                         avg_meth_sample = colMeans(ccle$meth, na.rm=TRUE))
cdf = ccle$clines %>%
    select(CCLE_ID, Name, Site_Primary, cohort=tcga_code) %>%
    filter(cohort %in% args$tissue) %>%
    left_join(avg_meth_sample) %>%
    left_join(reshape2::melt(ccle$copies, value.name="copies")) %>%
    left_join(reshape2::melt(ccle$meth, value.name="meth")) %>%
    filter(copies > 2-et) %>%
    group_by(gene) %>%
        tidyr::nest() %>%
    ungroup()

res = cdf %>%
    mutate(res = clustermq::Q(fit_beta, df=data, n_jobs=as.integer(args$cores),
        memory=as.integer(args$memory), pkgs=c("dplyr", "betareg"))) %>%
    tidyr::unnest("res") %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value) %>%
    select(-data) %>%
    filter(gene != "PSMB4") # p ~ 0, why?

pdf(args$plotfile, 10, 8)
print(plt$volcano(res, text.size=2.5, label_top=20, pos_label_bias=0.2))
dev.off()

writexl::write_xlsx(res, args$outfile)
