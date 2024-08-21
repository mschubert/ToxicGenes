library(dplyr)
library(ggplot2)
cm = import('../../report/common')

dset = readr::read_tsv("../../cor_tcga_ccle/positive_comp_set.tsv")

ours = c(na.omit(dset$gene[with(dset, est_ccle < -0.3 & est_tcga < -0.3)]))
goncalves = readxl::read_xlsx("mmc3.xlsx", 2) %>%
    filter(Transcriptomics_Proteomics > 0.3) %>%
    pull(Gene)
schukken = readxl::read_xlsx("Supplemental_Table_S2.xlsx", sheet="Difference upon aneuploidy")
sch_gene = schukken %>%
    filter(RNA.Gain.Difference.Category %in% c("Buffering", "Anti-Scaling")) %>%
    pull(Gene_Symbol)
sch_prot = schukken %>%
    filter(Protein.Gain.Difference.Category %in% c("Buffering", "Anti-Scaling")) %>%
    pull(Gene_Symbol)

ds = dset %>% select(gene, stat_orf) %>% na.omit()
ds = ds %>% mutate(type = "All genes") %>%
    bind_rows(ds %>% filter(gene %in% goncalves) %>% mutate(type="Goncalves")) %>%
    bind_rows(ds %>% filter(gene %in% sch_gene) %>% mutate(type="Schukken_Gene")) %>%
    bind_rows(ds %>% filter(gene %in% sch_prot) %>% mutate(type="Schukken_Protein")) %>%
    bind_rows(ds %>% filter(gene %in% ours) %>% mutate(type="Ours")) %>%
    mutate(type = factor(type, levels=unique(type)))

p = ggplot(ds, aes(x=type, y=stat_orf, fill=type)) +
    geom_boxplot(outlier.shape=NA, alpha=0.7) +
    ggsignif::geom_signif(y_position=c(3.2, 4.4, 5.5, 6.6), color="black", test=t.test,
        map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
        comparisons=list(c("All genes", "Goncalves"),
                         c("All genes", "Schukken_Gene"),
                         c("All genes", "Schukken_Protein"),
                         c("All genes", "Ours"))) +
    coord_cartesian(ylim=c(-7.5, 9)) +
    labs(fill = "Study", x = "Study", y = "Δ ORF (Wald statistic)") +
#    scale_fill_manual(values=cm$cols[c("Background", "Compensated", "Hyperactivated")]) +
    theme_classic() +
    theme(axis.text.x = element_blank()) +
    geom_hline(yintercept=median(ds$stat_orf[ds$type=="All genes"], na.rm=TRUE),
               linetype="dashed", color="black")

ll = list(ours=ours, Goncalves=goncalves, `Schukken gene`=sch_gene, `Schukken\nprotein`=sch_prot) |>
    lapply(function(x) na.omit(x[!duplicated(x)]))

saveRDS(list(overlap=ll, genes=ds), file="compensation.rds")

pdf("compensation.pdf", 5, 4)
print(p)
print(ggvenn::ggvenn(ll, set_name_size=4, text_size=3))
dev.off()
