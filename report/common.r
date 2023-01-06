import_package("dplyr", attach=TRUE)

hlg = c("MYC", "EGFR", "CCND1", "CDKN1A", "TP53", "BAP1", "CDKN1A", "IL7R", "CKS1B",
        "APC", "CDKN2A", "KRAS", "NRAS", "RB1", "CCNE1", "PIK3CA", "AURKA")

get_cosmic_annot = function() {
    manual = readRDS("../data/genesets/manual.rds")
    manual[grepl("Cosmic_(OG|TSG)_Tier", names(manual))] %>%
        stack() %>% as_tibble() %>%
        transmute(gene_name = values,
                  type = case_when(
                      grepl("OG", ind) ~ "Oncogene",
                      grepl("TSG", ind) ~ "TSG"
                  ),
                  tier = sub(".*Tier([12])$", "\\1", ind)) %>%
        distinct() %>%
        group_by(gene_name, tier) %>%
        summarize(type = ifelse(length(type) == 1, as.character(type), "OG+TSG"))
}

cols = c(
    Genes="grey", Background="grey", Euploid="grey",
    Amplification="#b06166", Deletion="#8484A8",
    Amplified="#b06166", Deleted="#8484A8", "Amp+Del"="#f0f0f0",
    Oncogene="#9f7ccd", TSG="#4ec472", "OG+TSG"="#d59e88",
    Compensated="#de493d", Hyperactivated="#adc5ee", Hyp2="#7aa1e3",
    CCLE="#226b94", TCGA="#74ad9b", ORF="#f7974e",
    "ORF dropout"="#f7974e",
    "Comp+ORF"="#E35740"
)

fmt_p = function(p) {
    lab = sprintf("italic(P)~`=`~%.2g", p)
    sub("[0-9.]+e", "10^", lab)
}
