library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'orf_naive', 'xlsx', '../orf/fits_naive.xlsx'),
    opt('z', 'orf_corr', 'xlsx', '../orf/fits_corrected.xlsx'),
    opt('c', 'ccle', 'xlsx', '../ccle/pan/stan-nb.xlsx'),
    opt('t', 'tcga_naive', 'xlsx', '../tcga/pan/stan-nb_naive.xlsx'),
    opt('p', 'tcga_pur', 'xlsx', '../tcga/pan/stan-nb_pur.xlsx'),
    opt('l', 'plotfile', 'pdf', 'pan.pdf')
)

orf1 = readxl::read_xlsx(args$orf_naive) %>% dplyr::rename(gene=`GENE SYMBOL`)
orf2 = readxl::read_xlsx(args$orf_corr) %>% dplyr::rename(gene=`GENE SYMBOL`)
ccle = readxl::read_xlsx(args$ccle) %>%
    filter(p.value < 1e-3) %>% mutate(estimate = pmax(-2, pmin(estimate, 2.5)))
tcga1 = readxl::read_xlsx(args$tcga_naive) %>%
    filter(p.value < 1e-5) %>% mutate(estimate = pmax(-2, pmin(estimate, 2.5)))
tcga2 = readxl::read_xlsx(args$tcga_pur) %>%
    filter(p.value < 1e-5) %>% mutate(estimate = pmax(-2, pmin(estimate, 2.5)))

hl = c("RBM12", "RBM14", "SNRPA", "HNRNPL", "CDKN1A", "H3F3C", "DAZAP1", "BANP", "HHEX",
       "CEBPE", "HNRNPL", "HNRNPA2B1", "KLF2", "KLF4", "IRF2", "SOX15", "ELK3", "ADORA2A",
       "BCL6", "KLF12", "IFNG", "IFNA", "ZBTB14", "FOS", "PIK3IP1", "RBMS2", "KLF4", "HES1",
       "LCOR", "ZC3H10", "ALKBH7", "CREB1", "GPR3", "KLF3")
hl2 = grep("^(TR[ABG]V|IG[KL]V)", tcga1$gene, value=T)

# tcga purity
both = inner_join(tcga1 %>% transmute(gene=gene, tcga_naive=estimate),
           tcga2 %>% transmute(gene=gene, tcga_pur=sign(estimate)*pmin(abs(estimate), 5)))
p1 = ggplot(both, aes(x=tcga_naive, y=tcga_pur)) +
    geom_point(alpha=0.1) +
    geom_point(data=both %>% filter(gene %in% hl), color="red") +
    geom_point(data=both %>% filter(gene %in% hl2), color="blue") +
    geom_text(data=both %>% filter(gene %in% hl2), aes(label=gene), size=1.5) +
    geom_density2d(color="green", size=0.2, breaks=c(0.05, 0.1,0.2,0.9)) +
    ggrepel::geom_label_repel(data=both %>% filter(gene %in% hl), aes(label=gene), size=2.5)

# tcga-ccle cor
both = ccle %>% transmute(gene=gene, ccle=estimate) %>%
    left_join(tcga1 %>% transmute(gene=gene, tcga_naive=estimate, tcga_aneup=n_aneup)) %>%
    left_join(tcga2 %>% transmute(gene=gene, tcga_pur=estimate))
m0 = broom::glance(lm(ccle ~ tcga_naive, data=both)); m0
p2 = ggplot(both, aes(x=ccle, y=tcga_naive)) +
    geom_point(aes(size=pmin(500, tcga_aneup)), alpha=0.1) +
    geom_point(data=both %>% filter(gene %in% hl), aes(size=pmin(500, tcga_aneup)), color="red") +
    scale_size_area() +
    geom_density2d(color="green", size=0.2, breaks=c(0.05, 0.1,0.2,0.9)) +
    geom_smooth(method="lm") +
    ggrepel::geom_label_repel(data=both %>% filter(gene %in% hl), aes(label=gene), size=2.5, max.overlaps=Inf) +
    labs(title = "stan-nb", subtitle=sprintf("p=%.2g (R^2=%.2g)", m0$p.value, m0$r.squared))

m0 = broom::glance(lm(ccle ~ tcga_pur, data=both)); m0
p3 = ggplot(both, aes(x=ccle, y=tcga_pur)) +
    geom_point(aes(size=pmin(500, tcga_aneup)), alpha=0.1) +
    geom_point(data=both %>% filter(gene %in% hl), aes(size=pmin(500, tcga_aneup)), color="red") +
    scale_size_area() +
    geom_density2d(color="green", size=0.2, breaks=c(0.05, 0.1,0.2,0.9)) +
    geom_smooth(method="lm") +
    ggrepel::geom_label_repel(data=both %>% filter(gene %in% hl), aes(label=gene), size=2.5, max.overlaps=Inf) +
    labs(title = "stan-nb", subtitle=sprintf("p=%.2g (R^2=%.2g)", m0$p.value, m0$r.squared))

# tcga+ccle vs orf screen
cmp = both %>% mutate(gene=gene, tcga_ccle = (tcga_pur + ccle)/2) %>%
#    inner_join(orf1 %>% transmute(gene=gene, orf1=estimate)) %>%
#    inner_join(orf2 %>% transmute(gene=gene, orf2=estimate))
    inner_join(orf1 %>% transmute(gene=gene, orf1=statistic)) %>%
    inner_join(orf2 %>% transmute(gene=gene, orf2=statistic))
m1 = broom::glance(lm(orf1 ~ tcga_ccle, data=cmp)); m1
m2 = broom::glance(lm(orf2 ~ tcga_ccle, data=cmp)); m2
p4 = ggplot(cmp, aes(x=tcga_ccle, y=orf2)) +
    geom_point(aes(size=pmin(500, tcga_aneup)), alpha=0.1) +
    scale_size_area() +
    geom_density2d(color="green", size=0.2, breaks=c(0.05, 0.1,0.2,0.9)) +
    geom_smooth(method="lm") +
    geom_point(data=cmp %>% filter(gene %in% hl), aes(size=pmin(500, tcga_aneup)), color="red") +
    ggrepel::geom_label_repel(data=cmp %>% filter(gene %in% hl), #filter(tcga_ccle < -0.5, orf2 < -1),
                              aes(label=gene), size=2.5, max.overlaps = Inf) +
    labs(title = "stan-nb", subtitle=sprintf("p=%.2g [R^2 = %.2g (z), %.2g (loess z)]",
                                             m1$p.value, m1$r.squared, m2$r.squared))

pdf(args$plotfile)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
