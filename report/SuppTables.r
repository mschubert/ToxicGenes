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

writexl::write_xlsx(ccle[c("Pan-Cancer", sort(setdiff(names(ccle), "Pan-Cancer")))], "TableS1_CCLE-comp.xlsx")
writexl::write_xlsx(tcga[c("Pan-Cancer", sort(setdiff(names(tcga), "Pan-Cancer")))], "TableS2_TCGA-comp.xlsx")
writexl::write_xlsx(orf, "TableS3_ORF-toxicity.xlsx")
