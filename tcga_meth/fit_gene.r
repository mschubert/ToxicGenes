library(dplyr)
library(betareg)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')

fit_beta = function(df) {
    df$meth = pmin(df$meth, 1 - 1e-3)
    df$meth = pmax(df$meth, 1e-3)

    if (length(unique(df$cohort)) == 1)
        fml = meth ~ purity + avg_meth_sample + copies
    else
        fml = meth ~ cohort + purity + avg_meth_sample + copies

    tryCatch({
            betareg(fml, data=df) %>%
                broom::tidy() %>%
                filter(term == "copies") %>%
                select(-component, -term)
        }, error = function(e) data.frame(estimate=NA))
}

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('t', 'tissue', 'pan|TCGA cohort', 'pan'),
    opt('e', 'euploid_tol', 'numeric copy dev from euploid', '0.15'),
    opt('c', 'cores', 'parallel cores', '5'),
    opt('m', 'memory', 'mem per core', '1024'),
    opt('o', 'outfile', 'rds', 'pan/betareg.rds'),
    opt('p', 'plotfile', 'pdf', 'pan/betareg.pdf')
)

cfg = yaml::read_yaml(args$config)
et = as.numeric(args$euploid_tol)
if (args$tissue == "pan")
    args$tissue = setdiff(cfg$tcga_tissues, "pan") # reduced set

copies = lapply(args$tissue, tcga$cna_genes, gene="external_gene_name") %>%
    narray::stack(along=2)
meth = tcga$meth(args$tissue, cpg="avg")
names(dimnames(copies)) = names(dimnames(meth)) = c("gene", "sample")
avg_meth_sample = tibble(sample = colnames(meth),
                         cohort = tcga$barcode2study(colnames(meth)),
                         avg_meth_sample = colMeans(meth, na.rm=TRUE))
purity = tcga$purity() %>% select(sample=Sample, cohort, purity=estimate) %>% na.omit()

cdf = purity %>%
    inner_join(avg_meth_sample) %>%
    left_join(reshape2::melt(copies, value.name="copies") %>% filter(copies > 2-et)) %>%
    left_join(reshape2::melt(meth, value.name="meth")) %>%
    na.omit() %>%
    group_by(gene) %>%
        tidyr::nest() %>%
    ungroup()

res = cdf %>%
    mutate(res = clustermq::Q(fit_beta, df=data, n_jobs=as.integer(args$cores),
        memory=as.integer(args$memory), pkgs=c("dplyr", "betareg"))) %>%
    select(-data) %>%
    tidyr::unnest("res") %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

pdf(args$plotfile, 10, 8)
print(plt$volcano(res, text.size=2.5, label_top=20, pos_label_bias=0.2))
dev.off()

writexl::write_xlsx(res, args$outfile)
