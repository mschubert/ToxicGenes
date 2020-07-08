library(dplyr)
library(ggplot2)
tcga = import('data/tcga')
util = import('./util')
#rlm3 = import('../tcga/fit_rlm3')

mut = tcga$mutations() %>%
    filter(Study == "BRCA", Hugo_Symbol == "TP53") %>%
    transmute(Sample = Tumor_Sample_Barcode,
              Variant = Variant_Classification,
              Var = Variant_Type,
              AAChange = AAChange)

snp = mut %>% filter(Var == "SNP") %>% pull(Sample) %>% unique()
del = mut %>% filter(Var == "DEL") %>% pull(Sample) %>% unique()

td = util$load_tcga("BRCA", top="CDKN1A") %>% #FIXME: is expr libsize corrected?
    select(-p53_mut) %>%
    mutate(p53_mut = case_when(
        sample %in% snp ~ "p53_snp",
        sample %in% del ~ "p53_del",
        TRUE ~ "p53_wt")) %>%
    filter(cancer_copies >= 2-0.15)
table(td$p53_mut)
#rlm3$do_fit()

mean_eup = td %>% filter(cancer_copies < 2+0.15) %>%
    group_by(gene, p53_mut) %>%
    summarize(mean_eup = mean(expr))

pdf("BRCA-p53-p53.pdf", 10, 4)
ggplot(td, aes(x=cancer_copies, y=expr)) +
    geom_abline(data=mean_eup, aes(slope=mean_eup/2, intercept=0), linetype="dashed", color="red") +
    geom_point(aes(color=meth_eup_scaled), size=2, alpha=0.5) +
    geom_smooth(method="lm") +
    scale_color_distiller(palette="RdBu") +
    facet_grid(gene ~ p53_mut)

ggplot(td, aes(x=cancer_copies, y=meth)) +
    geom_point(size=2) +
    geom_smooth(method="lm") +
    facet_grid(gene ~ p53_mut)
dev.off()
