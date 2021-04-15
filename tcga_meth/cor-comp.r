library(dplyr)
library(ggplot2)
sys = import('sys')

join_plot_one = function(comp) {
    both = comp %>%
        transmute(gene = gene, comp=estimate) %>%
        inner_join(meth) %>%
        left_join(orf)

    broom::tidy(lm(meth ~ comp, data=both))

    ggplot(data=both %>% filter(!is.na(orf_dir)), aes(x=comp, y=meth)) +
        geom_point(data=both %>% filter(is.na(orf_dir)), alpha=0.1) +
        geom_point(aes(color=orf_dir), alpha=0.1) +
        geom_density2d(aes(color=orf_dir)) +
        geom_smooth(method="lm")
}

args = sys$cmd$parse(
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('n', 'comp_naive', 'rds', '../tcga/pan/stan-nb_naive.rds'),
    opt('u', 'comp_pur', 'rds', '../tcga/pan/stan-nb_pur.rds'),
    opt('a', 'comp_puradj', 'rds', '../tcga/pan/stan-nb_puradj.rds'),
    opt('o', 'orf', 'xlsx', '../orf/pan/genes.xlsx'),
    opt('m', 'meth', 'xlsx', 'pan/gene.xlsx'),
    opt('p', 'plotfile', 'pdf', 'pan/cor-comp.pdf')
)

comp_naive = readRDS(args$comp_naive)
comp_pur = readRDS(args$comp_pur)
comp_puradj = readRDS(args$comp_puradj) %>%
    mutate(estimate = sign(estimate) * pmin(abs(estimate), 2))

meth = readxl::read_excel(args$meth) %>%
    transmute(gene=gene, meth=statistic)
orf = readxl::read_excel(args$orf) %>%
    filter(adj.p < 0.1) %>%
    transmute(gene=name, orf=statistic, orf_dir=factor(sign(statistic)))

pdf(args$plotfile)
join_plot_one(comp_naive) + ggtitle("naive")
join_plot_one(comp_pur) + ggtitle("pur")
join_plot_one(comp_puradj) + ggtitle("puradj")
dev.off()
