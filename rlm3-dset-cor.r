library(dplyr)
library(ggplot2)

orf1 = readxl::read_xlsx("orf/fits_naive.xlsx") %>% dplyr::rename(gene=`GENE SYMBOL`)
orf2 = readxl::read_xlsx("orf/fits_corrected.xlsx") %>% dplyr::rename(gene=`GENE SYMBOL`)
ccle = readxl::read_xlsx("ccle/pan/rlm3.xlsx", "amp") %>% dplyr::rename(gene=name) %>%
    filter(adj.p < 1) %>% mutate(estimate = pmax(-2, pmin(estimate, 2.5)))
tcga1 = readxl::read_xlsx("tcga/pan/rlm3_naive.xlsx", "amp") %>% dplyr::rename(gene=name) %>%
    filter(adj.p < 1) %>% mutate(estimate = pmax(-2, pmin(estimate, 2.5)))
tcga3 = readxl::read_xlsx("tcga/pan/rlm3_pur.xlsx", "amp") %>% dplyr::rename(gene=name) %>%
    filter(adj.p < 1) %>% mutate(estimate = pmax(-2, pmin(estimate, 2.5)))
#tcga3 = readxl::read_xlsx("tcga/pan/rlm3_puradj.xlsx", "amp")

hl = c("RBM12", "RBM14", "SNRPA", "HNRNPL", "CDKN1A", "H3F3C", "DAZAP1", "BANP", "HHEX",
       "CEBPE", "HNRNPL", "HNRNPA2B1", "KLF2", "KLF4", "IRF2", "SOX15", "ELK3", "ADORA2A",
       "BCL6", "KLF12", "IFNG", "IFNA", "ZBTB14", "FOS", "PIK3IP1", "RBMS2", "KLF4", "HES1",
       "LCOR", "ZC3H10", "ALKBH7", "CREB1", "GPR3", "KLF3")
hl = c(hl, grep("^(TR[ABG]V|IG[KL]V)", tcga1$gene, value=T))

# tcga purity
both = inner_join(tcga1 %>% transmute(gene=gene, tcga_naive=estimate),
           tcga3 %>% transmute(gene=gene, tcga_pur=estimate))
ggplot(both, aes(x=tcga_naive, y=tcga_pur)) +
    geom_point(alpha=0.1) +
    geom_point(data=both %>% filter(gene %in% hl), color="red") +
    ggrepel::geom_label_repel(data=both %>% filter(gene %in% hl), aes(label=gene), size=2.5)

# tcga-ccle cor
both = inner_join(ccle %>% transmute(gene=gene, ccle=estimate),
           tcga1 %>% transmute(gene=gene, tcga=estimate, tcga_aneup=n_aneup))
m0 = broom::glance(lm(ccle ~ tcga, data=both)); m0
ggplot(both, aes(x=ccle, y=tcga)) +
    geom_point(aes(size=pmin(500, tcga_aneup)), alpha=0.1) +
    geom_point(data=both %>% filter(gene %in% hl), aes(size=pmin(500, tcga_aneup)), color="red") +
    scale_size_area() +
    ggrepel::geom_label_repel(data=both %>% filter(gene %in% hl), aes(label=gene), size=2.5, max.overlaps=Inf) +
    labs(title = "rlm3", subtitle=sprintf("R^2=%.2g", m0$r.squared))

# tcga+ccle vs orf screen
cmp = both %>% mutate(gene=gene, tcga_ccle = (tcga + ccle)/2) %>%
    inner_join(orf1 %>% transmute(gene=gene, orf1=estimate)) %>%
    inner_join(orf2 %>% transmute(gene=gene, orf2=estimate))
m1 = broom::glance(lm(orf1 ~ tcga_ccle, data=cmp)); m1
m2 = broom::glance(lm(orf2 ~ tcga_ccle, data=cmp)); m2
ggplot(cmp, aes(x=tcga_ccle, y=orf2)) +
    geom_point(aes(size=pmin(500, tcga_aneup)), alpha=0.1) +
    scale_size_area() +
    geom_smooth(method="lm") +
    geom_point(data=cmp %>% filter(gene %in% hl), aes(size=pmin(500, tcga_aneup)), color="red") +
    ggrepel::geom_label_repel(data=cmp %>% filter(gene %in% hl), #filter(tcga_ccle < -0.5, orf2 < -1),
                              aes(label=gene), size=2.5, max.overlaps = Inf)
