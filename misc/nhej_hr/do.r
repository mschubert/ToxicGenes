library(dplyr)
library(ggplot2)
library(patchwork)
gset = import('genesets')

go = gset$get_human("GO_Biological_Process_2021")
hr = grep("homologous recomb", names(go), value=TRUE)
nh = grep("nonhom", names(go), value=TRUE)
gs = tibble(type = c(rep("HR", length(hr)), rep("NHEJ", length(nh))),
            name = c(hr, nh),
            gene = go[c(hr, nh)]) %>%
    tidyr::unnest(gene)
gs$name = factor(gs$name, levels=unique(gs$name))

dset = readr::read_tsv("../../cor_tcga_ccle/positive_comp_set.tsv") %>%
    mutate(is_orf_hit = stat_orf < -5 & est_orf < log2(0.7) & !is.na(stat_orf)) %>%
    filter(gene %in% gs$gene) %>%
    mutate(gene = forcats::fct_reorder(gene, est_tcga)) %>%
    select(gene, CCLE=est_ccle, TCGA=est_tcga, ORF=stat_orf) %>%
    tidyr::gather("dataset", "value", -gene) %>%
    mutate(dataset = factor(dataset, levels=c("TCGA", "CCLE", "ORF")))

gs$gene = factor(gs$gene, levels=levels(dset$gene))
gs = na.omit(gs)

p1 = ggplot(gs, aes(y=gene, x=name, fill=type)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=1)) +
    scale_fill_brewer(palette="Dark2")

p2 = ggplot(dset, aes(y=gene, x=value, fill=dataset)) +
    geom_col() +
    facet_grid(. ~ dataset, scales="free")

pdf("HR_vs_NHEJ.pdf", 12,25)
(p1 | p2) + plot_layout(widths=c(1,5), guides="collect")
dev.off()
