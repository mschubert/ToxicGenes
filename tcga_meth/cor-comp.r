library(dplyr)
library(ggplot2)
sys = import('sys')

join_plot_one = function(comp) {
    both = comp %>%
#        mutate(estimate = estimate * (1-adj.p)) %>%
        transmute(gene = gene, comp=sign(estimate) * pmin(abs(estimate), 2), n_aneup=n_aneup) %>%
        filter(n_aneup >= min_aneup) %>%
        inner_join(meth) %>%
        left_join(orf)

    m1 = broom::tidy(lm(meth ~ comp, data=both)) %>% filter(term == "comp")
    m2 = broom::tidy(lm(meth ~ comp, data=both %>% filter(orf_dir == 1))) %>% filter(term == "comp") # hyperact
    m3 = broom::tidy(lm(meth ~ comp, data=both %>% filter(orf_dir == -1))) %>% filter(term == "comp") # comp

    ggplot(data=both %>% filter(!is.na(orf_dir)), aes(x=comp, y=meth)) +
        geom_point(aes(size=n_aneup), data=both %>% filter(is.na(orf_dir)), alpha=0.05) +
        geom_point(aes(size=n_aneup, color=orf_dir), alpha=0.3) +
        geom_density2d(aes(color=orf_dir), size=0.1, show.legend=FALSE) +
        geom_density2d(color="black", size=0.1, show.legend=FALSE) +
        geom_smooth(method="lm", se=FALSE, color="black") +
        geom_smooth(aes(color=orf_dir), method="lm", se=FALSE) +
        ggrepel::geom_text_repel(aes(label=gene), data=both %>% filter(orf_dir == -1), size=1.8,
                                 max.overlaps=20, segment.alpha=0.3) +
        scale_size(range=c(0.1, 5)) +
        labs(subtitle = sprintf("%.2f (p=%.2g) [hyper: %.2f (p=%.2g) | comp: %.2f (p=%.2g)]",
                                m1$estimate, m1$p.value,
                                m2$estimate, m2$p.value,
                                m3$estimate, m3$p.value))
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

if (args$tissue == "pan") {
    min_aneup = 100
} else {
    min_aneup = 20
}

comp_naive = readRDS(args$comp_naive)
comp_pur = readRDS(args$comp_pur)
comp_puradj = readRDS(args$comp_puradj)

meth = readxl::read_excel(args$meth) %>%
    transmute(gene=gene, meth=estimate*(1-adj.p))
orf = readxl::read_excel(args$orf) %>%
    filter(adj.p < 0.1) %>%
    transmute(gene=name, orf=statistic, orf_dir=factor(sign(statistic)))

pdf(args$plotfile)
join_plot_one(comp_naive) + ggtitle("naive")
join_plot_one(comp_pur) + ggtitle("pur")
join_plot_one(comp_puradj) + ggtitle("puradj")
dev.off()
