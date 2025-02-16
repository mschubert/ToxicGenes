import_package("ggplot2", attach=TRUE)
import_package("dplyr", attach=TRUE)

hlg = c("MYC", "EGFR", "CCND1", "CDKN1A", "TP53", "BAP1", "CDKN1A", "IL7R", "CKS1B",
        "APC", "CDKN2A", "KRAS", "NRAS", "RB1", "CCNE1", "PIK3CA", "AURKA",
        "FGFR1", "NOTCH2")

get_cosmic_annot = function() {
    manual = readRDS(module_file("../data/genesets/manual.rds"))
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
    cfg = yaml::read_yaml(module_file("../config.yaml"))
    fname = module_file("SuppData2_CCLE-comp.xlsx")
    ccle = readxl::excel_sheets(fname) %>%
        setdiff("description") %>%
        sapply(readxl::read_xlsx, path=fname, simplify=FALSE) %>%
        bind_rows(.id="tissue") %>%
        mutate(src = "CCLE")
    fname = module_file("SuppData3_TCGA-comp.xlsx")
    tcga = readxl::excel_sheets(fname) %>%
        setdiff("description") %>%
        sapply(readxl::read_xlsx, path=fname, simplify=FALSE) %>%
        bind_rows(.id="tissue") %>%
        mutate(src = "TCGA")

    bind_rows(ccle, tcga) %>%
        group_by(tissue, gene) %>%
            filter(all(c("CCLE", "TCGA") %in% src)) %>%
            mutate(is_common = all(type == "Compensated")) %>%
        ungroup() %>%
        mutate(tissue = factor(tissue) %>% relevel("Pan-Cancer"))
}

get_tox = function() {
    fname = module_file("SuppData5_ORF-toxicity.xlsx")
    readxl::excel_sheets(fname) %>%
        setdiff("description") %>%
        sapply(readxl::read_xlsx, path=fname, simplify=FALSE)
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
    intersect(comp, tox$gene)
}

get_pancan_summary = function() {
    ccle = readxl::read_xlsx("SuppData2_CCLE-comp.xlsx", sheet="Pan-Cancer") |>
        transmute(gene, comp_ccle=compensation, type_ccle=type)
    tcga = readxl::read_xlsx("SuppData3_TCGA-comp.xlsx", sheet="Pan-Cancer") |>
        transmute(gene, comp_tcga=compensation, type_tcga=type)
    orf = readxl::read_xlsx("SuppData5_ORF-toxicity.xlsx", sheet="Pan-Cancer") |>
        transmute(gene, est_orf=estimate, stat_orf=statistic, is_tox=is_toxic)
    inner_join(ccle, tcga) |>
        mutate(type = ifelse(type_ccle == type_tcga, type_ccle, "Background"),
               is_comp = type == "Compensated") |>
        left_join(orf) |>
        mutate(is_argos = is_comp & is_tox)
}

cols = c(
    Genes="grey", Background="grey", Euploid="grey", Other="grey",
    Amplification="#b06166", Deletion="#8484A8",
    Amplified="#b06166", Deleted="#8484A8", "Amp+Del"="#f0f0f0",
    Oncogene="#9f7ccd", TSG="#4ec472", "OG+TSG"="#d59e88",
    Compensated="#de493d", `Comp.`="#de493d", Hyperactivated="#adc5ee", Hyp2="#7aa1e3",
    CCLE="#226b94", TCGA="#74ad9b", ORF="#f7974e",
    "ORF dropout"="#f7974e",
    "Comp+ORF"="#E35740"
)
col_study = c(
    `All genes` = "#c1c1c1", ours = "#7f7fff", Goncalves = "#ffff7f",
    `Schukken\n(gene)` = "#7fff7f", `Schukken\n(protein)` = "#ff7f7f",
    Other = "#c1c1c1", Comp. = "#7f7fff"
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
