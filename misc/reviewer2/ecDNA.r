library(dplyr)
library(EnsDb.Hsapiens.v75) # same versions as Kim/Verhaak Nat Gen paper

loc2gene = function(loc) {
    rng = GenomicRanges::GRanges(loc)
    ov = plyranges::join_overlap_intersect(genes, rng)
    ov$symbol
}

edb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
genes = genes(edb)

fname = "/home/m.schubert/mila/data/public_data/41588_2020_678_MOESM2_ESM.xlsx"
dset = readxl::read_xlsx(fname, sheet="Supplementary Table 1") |>
    dplyr::filter(grepl("^TCGA-", sample_barcode)) |>
    mutate(amplicon_intervals = stringr::str_split(amplicon_intervals, ",")) |>
    tidyr::unnest(amplicon_intervals) |>
    mutate(genes = lapply(amplicon_intervals, loc2gene))

freqs = dset |>
    select(sample_barcode, amplicon_classification, genes) |>
    tidyr::unnest(genes) |>
    group_by(amplicon_classification, genes) |>
    summarize(n = n_distinct(sample_barcode)) |>
    arrange(desc(n))
