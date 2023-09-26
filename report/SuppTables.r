library(dplyr)

load_comp = function(fname) {
    dset = readRDS(fname) %>%
        mutate(compensation = (1 - p.value) * estimate)
}

ccle = list.files("../model_compensation/fit_ccle-amp", recursive=TRUE, full.names=TRUE) %>%
    setNames(tools::file_path_sans_ext(basename(.))) %>% lapply(load_comp)
names(ccle)[names(ccle) == "pan"] = "Pan-Cancer"

tcga = list.files("../model_compensation/fit_tcga_puradj-amp", recursive=TRUE, full.names=TRUE) %>%
    setNames(tools::file_path_sans_ext(basename(.))) %>% lapply(load_comp)
names(tcga)[names(tcga) == "pan"] = "Pan-Cancer"

writexl::write_xlsx(ccle[c("Pan-Cancer", sort(setdiff(names(ccle), "Pan-Cancer")))], "TableS1_CCLE-comp.xlsx")
writexl::write_xlsx(tcga[c("Pan-Cancer", sort(setdiff(names(tcga), "Pan-Cancer")))], "TableS2_TCGA-comp.xlsx")

pan = readxl::read_xlsx("../model_orf/fits_naive.xlsx", sheet="pan")
cline = sapply(readxl::excel_sheets("../model_orf/fits_per_screen.xlsx"), readxl::read_xlsx,
               path="../model_orf/fits_per_screen.xlsx", simplify=FALSE)
orf = c(`Pan-Cancer`=list(pan), cline) %>%
    lapply(. %>% dplyr::rename(gene = `GENE SYMBOL`) %>% filter(gene != "LOC254896"))
writexl::write_xlsx(orf, "TableS3_ORF-toxicity.xlsx")
