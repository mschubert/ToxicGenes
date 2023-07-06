library(dplyr)
library(ggplot2)

dset = readRDS("dset.rds")
copies = dset$copies[c("CDKN1A", "RBM14"),] %>%
    reshape2::melt() %>% as_tibble() %>%
    select(Gene_Symbol=Var1, cline=Var2, copies=value)
prot = readxl::read_xlsx("Table_S2_Protein_Quant_Normalized.xlsx", sheet=2) %>%
    filter(Gene_Symbol %in% c("CDKN1A", "RBM14"))

pmat = prot[,c("Gene_Symbol", grep("_TenPx[0-9]+$", colnames(prot), value=TRUE))] %>%
    tidyr::pivot_longer(-Gene_Symbol, names_to="cline", values_to="protein") %>%
    mutate(protein = as.numeric(protein),
           cline = sub("_TenPx[0-9]+$", "", cline)) %>%
    inner_join(copies)

ggplot(pmat, aes(x=copies, y=protein)) +
    geom_hline(yintercept=0, color="grey", size=1, linetype="dashed") +
    geom_point(alpha=0.8) +
    geom_smooth(method="lm", se=FALSE) +
    theme_minimal() +
    facet_wrap(~ Gene_Symbol)
