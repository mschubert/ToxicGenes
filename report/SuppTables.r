library(dplyr)

ccle = list.files("../model_compensation/fit_ccle-amp", recursive=TRUE, full.names=TRUE) %>%
    setNames(tools::file_path_sans_ext(basename(.))) %>%
    lapply(readRDS) %>%
    lapply(. %>% mutate(compensation = (1 - p.value) * estimate))
writexl::write_xlsx(ccle, "TableS1_CCLE-comp.xlsx")

tcga = list.files("../model_compensation/fit_tcga_puradj-amp", recursive=TRUE, full.names=TRUE) %>%
    setNames(tools::file_path_sans_ext(basename(.))) %>%
    lapply(readRDS) %>%
    lapply(. %>% mutate(compensation = (1 - p.value) * estimate))
writexl::write_xlsx(tcga, "TableS2_TCGA-comp.xlsx")

orf = readxl::read_xlsx("../model_orf/fits_naive.xlsx")
writexl::write_xlsx(orf, "TableS3_ORF-toxicity.xlsx")
