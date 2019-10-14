library(dplyr)
library(ggplot2)
sys = import('sys')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('c', 'cna', 'aod or doa', 'doa'),
    opt('t', 'tissue', 'pan or TCGA identifier', 'OV'),
    opt('p', 'plotfile', 'pdf', 'heatmap_OV.pdf'))

events = readr::read_tsv("../data/gistic/TCGA.pancan12.130815.zigg_events.160923.txt")
if (args$tissue != "pan")
    events = events %>% filter(tcga$barcode2study(sample) == args$tissue)

arms = events %>%
    filter(arm_length %in% c(1,2),
           amplitude >= 0.1) %>%
    group_by(arm_length, chr) %>%
        filter(base_start == min(base_start) | base_end == max(base_end)) %>%
    ungroup()

cnas = events %>%
    filter(event_type %in% c("doa", "aod"), sample %in% arms$sample) %>%
    group_by(sample) %>%
        mutate(covered = sum(abs(base_end - base_start))) %>%
    ungroup() %>%
    mutate(chr = factor(chr, levels=gtools::mixedsort(unique(chr))),
           sample = forcats::fct_reorder(sample, covered))

p = ggplot(cnas) +
    geom_segment(aes(x=base_start, xend=base_end, y=sample, yend=sample)) +
    facet_grid(event_type ~ chr, space="free", scales="free") +
    ggtitle(args$tissue) +
    theme(axis.text.y = element_blank())

pdf(args$plotfile, 18, 14)
print(p)
dev.off()
