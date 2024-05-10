library(dplyr)
library(ggplot2)
tcga = import('data/tcga')
cm = import('../../report/common')

gistic = lapply(tcga$cohorts(), tcga$cna_gistic, thresh=TRUE) |>
    narray::stack(along=2) |>
    reshape2::melt() |>
    select(gene=Var1, Sample=Var2, gistic=value) |>
    mutate(gistic = c("loss", "loss", "euploid", "gain", "gain")[gistic+3]) |>
    filter(gistic != "loss")

comp = cm$get_comp_genes(pan=TRUE)
argos = cm$get_argos(pan=TRUE)

muts = lapply(tcga$cohorts(), . %>% tcga$mutations() %>% select(Sample, gene=Hugo_Symbol, Consequence)) |>
    bind_rows() |>
    inner_join(gistic) |>
    filter(Consequence %in% c("3_prime_UTR_variant", "5_prime_UTR_variant", "frameshift_variant",
                              "intron_variant", "missense_variant", "stop_gained", "synonymous_variant"))
n_smp = n_distinct(muts$Sample)

freqs = muts |>
    group_by(gene, gistic, Consequence) |>
        summarize(freq = n_distinct(Sample) / n_smp) |>
    ungroup() |>
    mutate(class = case_when(gene %in% argos ~ "ARGOS",
                             gene %in% comp ~ "Compensated",
                             TRUE ~ "Other")) |>
    mutate(class = factor(class, levels=c("Other", "Compensated", "ARGOS")))

pdf("mut_enrich_types.pdf", 6, 18)
ggplot(freqs, aes(x=class, y=freq)) +
    geom_violin(outlier.shape=NA) +
    ggbeeswarm::geom_quasirandom(aes(shape=class, fill=class), size=2, alpha=0.8) +
    scale_shape_manual(values=c(Other=NA, Compensated=21, ARGOS=21)) +
    scale_fill_manual(values=c(Other=NA, Compensated="#74ad9b", ARGOS="#de493d")) +
    scale_y_log10() +
    labs(x = "Gene class",
         y = "TCGA mutation frequency") +
    ggsignif::geom_signif(color="black", y_position=c(-0.9,-1), test=t.test,
        map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
        comparisons = list(c("Other", "Compensated"), c("Compensated", "ARGOS"))) +
    facet_grid(Consequence ~ gistic)
dev.off()
