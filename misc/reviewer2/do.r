library(dplyr)
library(EnsDb.Hsapiens.v75) # same versions as Kim/Verhaak Nat Gen paper



fname = "/home/m.schubert/mila/data/public_data/41588_2020_678_MOESM2_ESM.xlsx"
dset = readxl::read_xlsx(fname, sheet="Supplementary Table 1") |>
    filter(grepl("^TCGA-", sample_barcode)) |>
    mutate(amplicon_intervals = stringr::str_split(amplicon_intervals, ",")) |>
    tidyr::unnest(amplicon_intervals)
