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

pan = readxl::read_xlsx("../model_orf/fits_naive.xlsx", sheet="pan")
cline = sapply(readxl::excel_sheets("../model_orf/fits_per_screen.xlsx"), readxl::read_xlsx,
               path="../model_orf/fits_per_screen.xlsx", simplify=FALSE)
orf = c(pan=list(pan), cline) %>%
    lapply(. %>% dplyr::rename(gene_name = `GENE SYMBOL`) %>% filter(gene_name != "LOC254896"))
writexl::write_xlsx(orf, "TableS3_ORF-toxicity.xlsx")
