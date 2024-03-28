library(dplyr)
library(ggplot2)
tcga = import('data/tcga')
cm = import('../../report/common')

amp_freq = readRDS("../../data/gistic_smooth.rds")$genes |>
    filter(type == "amplification") |>
    mutate(freq_amp = ifelse(frac >= 0.15, "Freq. amp.", "Not")) |>
    select(gene = gene_name, freq_amp)

comp = cm$get_comp_genes(pan=TRUE)
argos = cm$get_argos(pan=TRUE)

muts = lapply(tcga$cohorts(), . %>% tcga$mutations() %>% select(Sample, gene=Hugo_Symbol)) |>
    bind_rows() |>
    inner_join(amp_freq)
n_smp = n_distinct(muts$Sample)

freqs = muts |>
    group_by(gene, freq_amp) |>
        summarize(freq = n_distinct(Sample) / n_smp) |>
    ungroup() |>
    mutate(class = case_when(gene %in% argos ~ "ARGOS",
                             gene %in% comp ~ "Compensated",
                             TRUE ~ "Other")) |>
    mutate(class = factor(class, levels=c("Other", "Compensated", "ARGOS")))

pdf("mut_enrich.pdf", 8, 4)
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
    facet_wrap(~freq_amp)
dev.off()
