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

td = util$load_tcga("BRCA", top="CDKN1A") %>%
    select(-p53_mut) %>%
    mutate(p53_mut = case_when(
        sample %in% snp ~ "snp",
        sample %in% del ~ "del",
        TRUE ~ "wt")) %>%
    filter(cancer_copies >= 2-0.15)
table(td$p53_mut)
#rlm3$do_fit()

pdf("BRCA-p53-p53.pdf", 10, 4)
ggplot(td, aes(x=cancer_copies, y=expr)) +
    geom_point(alpha=0.2) +
    geom_smooth(method="lm") +
    facet_grid(gene ~ p53_mut)
dev.off()
