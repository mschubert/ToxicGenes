library(dplyr)
library(ggplot2)
cm = import('../../report/common')

# this should be possible in one go instead of lapply
loc2gene = function(loc) {
    rng = GenomicRanges::GRanges(loc)
    ov = plyranges::join_overlap_intersect(genes, rng)
    ov$symbol
}

comp = cm$get_comp_tissue() |> filter(tissue == "Pan-Cancer") |> select(-tissue)

# same versions as Kim/Verhaak Nat Gen paper
edb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
genes = GenomicFeatures::genes(edb)

fname = "/home/m.schubert/mila/data/public_data/41588_2020_678_MOESM2_ESM.xlsx"
kim2020 = readxl::read_xlsx(fname, sheet="Supplementary Table 1") |>
    filter(grepl("^TCGA-", sample_barcode)) |>
    mutate(amplicon_intervals = stringr::str_split(amplicon_intervals, ",")) |>
    tidyr::unnest(amplicon_intervals) |>
    mutate(gene = lapply(amplicon_intervals, loc2gene))

n_smp = length(unique(kim2020$sample_barcode))
freqs = kim2020 |>
    select(sample_barcode, amplicon_classification, gene) |>
    tidyr::unnest(gene) |>
    filter(gene %in% comp$gene) |>
    group_by(amplicon_classification, gene) |>
        summarize(freq = n_distinct(sample_barcode) / n_smp) |>
        slice_max(n=100, order_by=freq) |>
    ungroup()

#dset = left_join(comp, freqs, relationship="many-to-many") |>
#    mutate(amplicon_classification = tidyr::replace_na(amplicon_classification, "All other"))
dset = inner_join(comp, freqs, relationship="many-to-many")
saveRDS(dset, file="ecDNA.rds")

pdf("ecDNA.pdf", 8, 4)
ggplot(dset, aes(x=amplicon_classification, y=compensation)) +
    geom_boxplot(outlier.shape=NA) +
    ggbeeswarm::geom_quasirandom(aes(fill=amplicon_classification), shape=21, alpha=0.5) +
    guides(fill="none") +
    labs(x = "Amplification class",
         y = "Compensation score") +
    facet_wrap(~ src) +
    ggsignif::geom_signif(y_position=c(0.7,0.85,1), color="black", test=wilcox.test,
        map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
        comparisons = list(c("BFB", "Circular"),
                           c("BFB", "Heavily-rearranged"),
                           c("BFB", "Linear")))
dev.off()
