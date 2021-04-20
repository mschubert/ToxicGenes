library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('c', 'comp', 'rds', '../ccle/pan/stan-nb/genes.rds'),
    opt('o', 'orf', 'xlsx', '../orf/pan/genes.xlsx'),
    opt('m', 'meth', 'xlsx', 'pan/gene.xlsx'),
    opt('p', 'plotfile', 'pdf', 'pan/cor-comp.pdf')
)

if (args$tissue == "pan") {
    min_aneup = 50
} else {
    min_aneup = 10
}

comp = readRDS(args$comp) %>%
    transmute(gene = gene, comp=sign(estimate) * pmin(abs(estimate), 2))
meth = readxl::read_excel("./pan/cpg.xlsx") %>%
    transmute(gene=gene, meth=statistic)
orf = readxl::read_excel(args$orf) %>%
    filter(adj.p < 0.1) %>%
    transmute(gene=name, orf=statistic, orf_dir=factor(sign(statistic)))

both = inner_join(comp, meth) %>%
    left_join(orf)

m1 = broom::tidy(lm(meth ~ comp, data=both)) %>% filter(term == "comp")
m2 = broom::tidy(lm(meth ~ comp, data=both %>% filter(orf_dir == 1))) %>% filter(term == "comp") # hyperact
m3 = broom::tidy(lm(meth ~ comp, data=both %>% filter(orf_dir == -1))) %>% filter(term == "comp") # comp

p = ggplot(data=both %>% filter(!is.na(orf_dir)), aes(x=comp, y=meth)) +
    geom_point(data=both %>% filter(is.na(orf_dir)), alpha=0.1) +
    geom_point(aes(color=orf_dir), alpha=0.1) +
    geom_density2d(aes(color=orf_dir)) +
    geom_smooth(method="lm") +
    geom_smooth(aes(color=orf_dir), method="lm", se=FALSE) +
    labs(subtitle = sprintf("%.2f (p=%.2g) [hyper: %.2f (p=%.2g) | comp: %.2f (p=%.2g)]",
                            m1$estimate, m1$p.value,
                            m2$estimate, m2$p.value,
                            m3$estimate, m3$p.value))

pdf(args$plotfile)
print(p)
dev.off()
