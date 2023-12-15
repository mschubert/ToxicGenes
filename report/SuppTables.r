library(dplyr)

load_comp = function(fname) {
    dset = readRDS(fname) %>%
        mutate(compensation = (1 - p.value) * estimate,
               is_comp = compensation < -0.3)
}

proc_tox = function(dset) {
    dset %>%
        dplyr::rename(gene = `GENE SYMBOL`) %>%
        filter(gene != "LOC254896") %>%
        mutate(is_toxic = p.value < 1e-5 & estimate < log2(0.7))
}

make_desc_comp = function(name) {
    list(description=tibble::tribble(
        name, ~ "",
        "", "",
        "gene", "HGNC gene symbol",
        "estimate", "The estimated compensation (beta 2) parameter without shrinkage",
        "std.error", "The standard error of the compensation estimate",
        "z_comp", "The z-score of the compensation estimate over the origin",
        "n_aneup", "Number of samples where gene is gained",
        "eup_reads", "Normalized rounded read count for all euploid samples",
        "n_eff", "MCMC effective sample size",
        "Rhat", "Markov chain convergence measure",
        "p.value", "The p-value calculated from the z-score",
        "adj.p", "The FDR-adjusted p-value",
        "compensation", "The compensation score (negative values mean compensated)",
        "is_comp", "Boolean value whether the gene is considered compensated",
    ))
}

make_desc_orf = function() {
    list(description=tibble::tribble(
        ~ "Supplementary Table 3: ORF Toxicity", ~ "",
        "", "",
        "gene", "HGNC gene symbol",
        "estimate", "The estimate of the linear regression model",
        "std.error", "The standard error of the estimate",
        "statistic", "Wald statistic",
        "p.value", "P-value of the regression model",
        "size", "How many data points were considered (cell lines x barcodes)",
        "adj.p", "FDR-adjusted p-value",
        "is_toxic", "Boolean value whether the gene is considered toxic"
    ))
}

ccle = list.files("../model_compensation/fit_ccle-amp", recursive=TRUE, full.names=TRUE) %>%
    setNames(tools::file_path_sans_ext(basename(.))) %>% lapply(load_comp)
names(ccle)[names(ccle) == "pan"] = "Pan-Cancer"

tcga = list.files("../model_compensation/fit_tcga_puradj-amp", recursive=TRUE, full.names=TRUE) %>%
    setNames(tools::file_path_sans_ext(basename(.))) %>% lapply(load_comp)
names(tcga)[names(tcga) == "pan"] = "Pan-Cancer"

pan = readxl::read_xlsx("../model_orf/fits_naive.xlsx", sheet="pan")
cline = sapply(readxl::excel_sheets("../model_orf/fits_per_screen.xlsx"), readxl::read_xlsx,
               path="../model_orf/fits_per_screen.xlsx", simplify=FALSE)
orf = c(`Pan-Cancer`=list(pan), cline) %>% lapply(proc_tox)
#tox_genes = with(orf$`Pan-Cancer`, gene[is_toxic])

#common = bind_rows(ccle, tcga, .id="src") %>%
#    select(gene, tissue, src) %>%
#    group_by(gene, tissue) %>%
#        mutate(is_common = all(is_comp),
#               is_argos = all(is_comp) & gene %in% tox_genes) %>%
#    ungroup()
#tcga = left_join(tcga, common)
#ccle = left_join(ccle, common)

writexl::write_xlsx(c(make_desc_comp(~ "Supplementary Table 1: CCLE compensation"),
                      ccle[c("Pan-Cancer", sort(setdiff(names(ccle), "Pan-Cancer")))]),
                    "TableS1_CCLE-comp.xlsx")
writexl::write_xlsx(c(make_desc_comp(~ "Supplementary Table 2: TCGA compensation"),
                      tcga[c("Pan-Cancer", sort(setdiff(names(tcga), "Pan-Cancer")))]),
                    "TableS2_TCGA-comp.xlsx")
writexl::write_xlsx(c(make_desc_orf(), orf), "TableS3_ORF-toxicity.xlsx")
