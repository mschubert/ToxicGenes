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

lookup = c(
    NCIH838_LUNG = "NCI-H838",
    NCIH1650_LUNG = "NCI-H1650",
    ZR751_BREAST = "ZR-75-1",
    HCC70_BREAST = "HCC70"
)
sel = pmat %>% filter(Gene_Symbol %in% c("CDKN1A", "RBM14")) %>%
    mutate(label = lookup[cline], has_label = !is.na(label))
cur = list(`Pan-can` = sel,
           Breast = sel %>% filter(grepl("BREAST", cline)),
           Lung = sel %>% filter(grepl("LUNG", cline))) %>%
    bind_rows(.id="Tissue") %>%
    mutate(Tissue = factor(Tissue, levels=c("Pan-can", "Breast", "Lung")))

plt = function(cur) {
    ggplot(cur, aes(x=pmin(copies, 5), y=pmin(2^protein, 5))) +
        geom_hline(yintercept=1, color="grey", size=1, linetype="dashed") +
        geom_point(color="grey", alpha=0.6) +
        geom_point(data=cur %>% filter(has_label), color="black", alpha=0.6) +
        geom_abline(slope=0.5, intercept=0, color="red", linetype="dashed", size=1) +
        geom_smooth(method="lm", se=FALSE, size=1) +
        ggrepel::geom_text_repel(aes(label=label)) +
        theme_minimal() +
        facet_grid(Gene_Symbol ~ Tissue) +
        labs(x = "copies",
             y = "Normalized protein expression")
}
pdf("protein_quant.pdf", 7.5, 3.5)
plt(cur %>% filter(Gene_Symbol == "CDKN1A"))
plt(cur %>% filter(Gene_Symbol == "RBM14"))
dev.off()
