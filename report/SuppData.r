library(dplyr)

load_comp = function(fname) {
    dset = readRDS(fname) %>%
        mutate(compensation = (1 - p.value) * estimate,
            type = case_when(
                compensation < -0.3 & adj.p <= 0.2 ~ "Compensated",
                compensation > 0.3 & adj.p <= 0.2 ~ "Hyperactivated",
                TRUE ~ NA_character_
            )
        ) %>%
        mutate(across(where(is.numeric) & !matches("p.value|adj.p"), ~ round(.x, 3)),
               across(matches("p.value|adj.p"), ~ sprintf("%.3g", .)))
}

proc_tox = function(dset) {
    dset %>%
        dplyr::rename(gene = `GENE SYMBOL`) %>%
        filter(gene != "LOC254896") %>%
        mutate(is_toxic = p.value < 1e-5 & estimate < log2(0.7)) %>%
        mutate(across(where(is.numeric) & !matches("p.value|adj.p"), ~ round(.x, 3)),
               across(matches("p.value|adj.p"), ~ sprintf("%.3g", .)))
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
        ~ "Supplementary Data 4: ORF Toxicity", ~ "",
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
    setNames(tools::file_path_sans_ext(basename(.))) %>% lapply(load_comp) %>%
    lapply(. %>% mutate(eup_reads=round(eup_reads), stroma_reads=round(stroma_reads)))
names(tcga)[names(tcga) == "pan"] = "Pan-Cancer"
# below: recover std.error from previous run where we haven't saved it, rerun is comp. expensive
tcga$`Pan-Cancer` = tcga$`Pan-Cancer` |> select(-n_genes) |> mutate(std.error = round(estimate/z_comp,3)) |>
    select(gene, estimate, std.error, everything())

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

fnames = c("ORF_DMSO_2019-02.txt", "ORF_DMSO-ETP_2019-07.xlsx", "20101003_ORF-size-plasmid.txt")
fpath = file.path("../data/orf/", fnames)
orfdata = list(readr::read_tsv(fpath[1]), readxl::read_xlsx(fpath[2]), readr::read_tsv(fpath[3])) |>
    setNames(tools::file_path_sans_ext(fnames))
libinfo = orfdata[[3]][9:10] |> filter(!(is.na(...9) & is.na(...10))) |>
    rename(lib=...9, cells=...10) |>
    tidyr::fill(lib, cells, .direction="down") |> filter(!is.na(cells))
libinfo$cells[9] = "ORF-T47D" # remove invalid character
orfdata$`Libraries` = libinfo
orfdata[[3]] = orfdata[[3]][1:7]

writexl::write_xlsx(c(make_desc_comp(~ "Supplementary Data 1: CCLE compensation"),
                      ccle[c("Pan-Cancer", sort(setdiff(names(ccle), "Pan-Cancer")))]),
                    "SuppData1_CCLE-comp.xlsx")
writexl::write_xlsx(c(make_desc_comp(~ "Supplementary Data 2: TCGA compensation"),
                      tcga[c("Pan-Cancer", sort(setdiff(names(tcga), "Pan-Cancer")))]),
                    "SuppData2_TCGA-comp.xlsx")
writexl::write_xlsx(orfdata, "SuppData3_ORFscreens.xlsx")
writexl::write_xlsx(c(make_desc_orf(), orf), "SuppData4_ORF-toxicity.xlsx")
