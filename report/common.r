import_package("ggplot2", attach=TRUE)
import_package("dplyr", attach=TRUE)

hlg = c("MYC", "EGFR", "CCND1", "CDKN1A", "TP53", "BAP1", "CDKN1A", "IL7R", "CKS1B",
        "APC", "CDKN2A", "KRAS", "NRAS", "RB1", "CCNE1", "PIK3CA", "AURKA",
        "FGFR1", "NOTCH2")

get_cosmic_annot = function() {
    manual = readRDS("../data/genesets/manual.rds")
    hms = manual[grepl("Cosmic_(OG|TSG)_Hallmark", names(manual))] %>%
        stack() %>% as_tibble() %>%
        transmute(gene_name = values,
                  type = case_when(
                      grepl("OG", ind) ~ "Oncogene",
                      grepl("TSG", ind) ~ "TSG"),
                  hallmark = TRUE)
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
            summarize(type = ifelse(length(unique(type)) == 1, as.character(type), "OG+TSG")) %>%
        ungroup() %>%
        left_join(hms) %>%
        mutate(hallmark = ifelse(is.na(hallmark), "No", "Yes"))
}

get_comp_tissue = function() {
    cfg = yaml::read_yaml("../config.yaml")
    ccle = readxl::excel_sheets("TableS1_CCLE-comp.xlsx") %>%
        sapply(readxl::read_xlsx, path="TableS1_CCLE-comp.xlsx", simplify=FALSE) %>%
        bind_rows(.id="tissue") %>%
        mutate(src = "CCLE")
    tcga = readxl::excel_sheets("TableS2_TCGA-comp.xlsx") %>%
        sapply(readxl::read_xlsx, path="TableS2_TCGA-comp.xlsx", simplify=FALSE) %>%
        bind_rows(.id="tissue") %>%
        mutate(src = "TCGA")

    bind_rows(ccle, tcga) %>%
        group_by(tissue, gene) %>%
            filter(all(c("CCLE", "TCGA") %in% src)) %>%
            mutate(is_common = all(is_comp)) %>%
        ungroup() %>%
        mutate(tissue = factor(tissue) %>% relevel("Pan-Cancer"))
}

get_tox = function() {
    sapply(readxl::excel_sheets("TableS3_ORF-toxicity.xlsx"), readxl::read_xlsx,
        path="TableS3_ORF-toxicity.xlsx", simplify=FALSE)
}

get_comp_genes = function(pan=FALSE) {
    comp = get_comp_tissue() %>% filter(is_common)
    if (pan)
        comp = comp[comp$tissue == "Pan-Cancer",]
    with(comp, intersect(gene[src=="CCLE"], gene[src=="TCGA"]))
}

get_argos = function(pan=FALSE) {
    comp = get_comp_genes(pan)
    tox = get_tox()$`Pan-Cancer` %>% filter(is_toxic)
    intersect(comp2, tox$gene)
}

cols = c(
    Genes="grey", Background="grey", Euploid="grey", Other="grey",
    Amplification="#b06166", Deletion="#8484A8",
    Amplified="#b06166", Deleted="#8484A8", "Amp+Del"="#f0f0f0",
    Oncogene="#9f7ccd", TSG="#4ec472", "OG+TSG"="#d59e88",
    Compensated="#de493d", Hyperactivated="#adc5ee", Hyp2="#7aa1e3",
    CCLE="#226b94", TCGA="#74ad9b", ORF="#f7974e",
    "ORF dropout"="#f7974e",
    "Comp+ORF"="#E35740"
)

fmt_p = function(p, sig=2) {
    lab = sprintf(paste0("italic(P)~`=`~%.", sig, "g"), p)
    sub("[0-9.]+e", "10^", lab)
}

text_sizes = function() {
    theme(plot.tag = element_text(size=24, face="bold"),
          plot.title = element_text(size=14),
          plot.subtitle = element_text(size=13),
          legend.title = element_text(size=12),
          legend.text = element_text(size=11),
          axis.title = element_text(size=12),
          axis.text = element_text(size=11),
          strip.text = element_text(size=12),
          strip.text.y.left = element_text(size=12),
          text = element_text(size=11),
          legend.background = element_blank(),
          legend.box.background = element_blank(),
          plot.background = element_blank())
}

theme_classic = function() ggplot2::theme_classic() + text_sizes()
theme_minimal = function() ggplot2::theme_minimal() + text_sizes()
