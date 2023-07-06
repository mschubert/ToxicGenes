library(dplyr)
library(ggplot2)

dset = readRDS("dset.rds")
copies = dset$copies %>%
    reshape2::melt() %>% as_tibble() %>%
    select(Gene_Symbol=Var1, cline=Var2, copies=value)
prot = readxl::read_xlsx("Table_S2_Protein_Quant_Normalized.xlsx", sheet=2)

pmat = prot[,c("Gene_Symbol", grep("_TenPx[0-9]+$", colnames(prot), value=TRUE))] %>%
    tidyr::pivot_longer(-Gene_Symbol, names_to="cline", values_to="protein") %>%
    mutate(protein = as.numeric(protein),
           cline = sub("_TenPx[0-9]+$", "", cline)) %>%
    inner_join(copies)

saveRDS(pmat, file="protein_quant.rds")

sel = pmat %>% filter(Gene_Symbol %in% c("CDKN1A", "RBM14"))
cur = list(`Pan-can` = sel,
           Breast = sel %>% filter(grepl("BREAST", cline)),
           Lung = sel %>% filter(grepl("LUNG", cline))) %>%
    bind_rows(.id="Tissue") %>%
    mutate(Tissue = factor(Tissue, levels=c("Pan-can", "Breast", "Lung")))

pdf("protein_quant.pdf", 8, 5)
ggplot(cur, aes(x=copies, y=protein)) +
    geom_hline(yintercept=0, color="grey", size=1, linetype="dashed") +
    geom_point(alpha=0.6) +
    geom_smooth(method="lm", se=FALSE) +
    theme_minimal() +
    facet_grid(Gene_Symbol ~ Tissue) +
    ylab("normalized protein expression")
dev.off()
