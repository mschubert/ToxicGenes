library(dplyr)

ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
    mutate(compensation = (1 - p.value) * estimate) %>%
    select(-n_genes)

tcga = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
    mutate(compensation = (1 - p.value) * estimate) %>%
    select(-n_genes)

orf = readxl::read_xlsx("../orf/fits_naive.xlsx")

dset = list(
    `S1 - CCLE compensation` = ccle,
    `S2 - TCGA compensation` = tcga,
    `S3 - ORF toxicity` = orf
)

writexl::write_xlsx(dset, "SuppTables.xlsx")
