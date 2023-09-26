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
    ccle = paste0(file.path("../model_compensation/fit_ccle-amp", cfg$ccle_tissues), ".rds") %>%
        lapply(readRDS) %>%
        setNames(cfg$ccle_tissues) %>% bind_rows(.id="tissue") %>%
        mutate(src = "CCLE")
    tcga = paste0(file.path("../model/fit_tcga_puradj-amp", cfg$tcga_tissues), ".rds") %>%
        lapply(readRDS) %>%
        setNames(cfg$tcga_tissues) %>% bind_rows(.id="tissue") %>%
        mutate(src = "TCGA")

    dset = bind_rows(ccle, tcga) %>%
        group_by(gene) %>% filter(all(c("CCLE", "TCGA") %in% src)) %>% ungroup() %>%
        mutate(tissue = ifelse(tissue == "pan", "Pan-Cancer", tissue),
               tissue = factor(tissue) %>% relevel("Pan-Cancer"),
               shrunk = estimate * (1-p.value))
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
