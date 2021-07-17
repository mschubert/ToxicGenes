library(dplyr)
library(ggplot2)
library(patchwork)
gdsc = import('data/gdsc')
tcga = import('data/tcga')

expr = t(gdsc$basal_expression()[c("FOXA1","MAX"),])
df = data.frame(cosmic_id = rownames(expr), expr) %>%
    as_tibble() %>%
    mutate(tissue = gdsc$cosmic$id2tissue(cosmic_id),
           cell_line = gdsc$cosmic$id2name(cosmic_id)) %>%
    filter(tissue %in% c("BRCA", "SKCM")) %>%
    mutate(label = ifelse(cell_line %in% c("MDA-MB-231", "HCC2157", "HCC1187",
        "MCF7", "HCC2218", "SKBR3", "MDA-MB-453", "T47D", "A2058", "COLO792",
        "HCC1954", "HCC1419", "MDA−MB−157", "BT−474"), cell_line, NA))

df2 = tidyr::gather(df, "gene", "expr", -cosmic_id, -tissue, -cell_line, -label)

p1 = ggplot(df2, aes(x=tissue, y=expr)) +
    geom_violin() +
    ggrepel::geom_text_repel(aes(label=label), size=3) +
    facet_wrap(~ gene) +
    ggtitle("GDSC")

expr2 = rbind(t(tcga$rna_seq("SKCM", trans="vst")[c("ENSG00000129514", "ENSG00000125952"),]),
              t(tcga$rna_seq("BRCA", trans="vst")[c("ENSG00000129514", "ENSG00000125952"),]))
colnames(expr2) = c("FOXA1", "MAX")
df3 = data.frame(cohort = tcga$barcode2study(rownames(expr2)), FOXA1=expr2[,"FOXA1"], MAX=expr2[,"MAX"]) %>%
    tidyr::gather("gene", "expr", -cohort)

p2 = ggplot(df3, aes(x=cohort, y=expr)) +
    geom_violin() +
    facet_wrap(~ gene) +
    ggtitle("TCGA")

p1 / p2
