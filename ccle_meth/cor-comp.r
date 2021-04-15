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

comp = readRDS(args$comp) %>%
    transmute(gene = gene, comp=estimate)
meth = readxl::read_excel("./pan/cpg.xlsx") %>%
    transmute(gene=gene, meth=statistic)
orf = readxl::read_excel(args$orf) %>%
    filter(adj.p < 0.1) %>%
    transmute(gene=name, orf=statistic, orf_dir=factor(sign(statistic)))

both = inner_join(comp, meth) %>%
    left_join(orf)

broom::tidy(lm(meth ~ comp, data=both))

pdf(args$plotfile)
ggplot(data=both %>% filter(!is.na(orf_dir)), aes(x=comp, y=meth)) +
    geom_point(data=both %>% filter(is.na(orf_dir)), alpha=0.1) +
    geom_point(aes(color=orf_dir), alpha=0.1) +
    geom_density2d(aes(color=orf_dir)) +
    geom_smooth(method="lm")
dev.off()
