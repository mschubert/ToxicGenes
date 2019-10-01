library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(plyranges)
seq = import('seq')
sys = import('sys')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('a', 'annot', 'rds w regions to check', 'exons.rds'),
    opt('n', 'num', 'number of top comp to plot', '20'),
    opt('c', 'cna', 'aod or doa', 'doa'),
    opt('t', 'tissue', 'pan or TCGA identifier', 'pan'),
    opt('p', 'plotfile', 'pdf', 'pan_doa.pdf'))

exons = readRDS(args$annot)

events = readr::read_tsv("../data/gistic/TCGA.pancan12.130815.zigg_events.160923.txt")
if (args$tissue != "pan")
    events = events %>% filter(tcga$barcode2study(sample) == "OV")
probes = readr::read_tsv("../data/gistic/genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1 2.txt",
                         col_names=FALSE) %>%
    dplyr::rename(id=X1, chr=X2, start=X3) %>%
    filter(! is.na(chr)) %>%
    makeGRangesFromDataFrame(end.field="start")

arms = events %>%
    filter(arm_length %in% c(1,2),
           amplitude >= 0.1) %>%
    group_by(arm_length, chr) %>%
        filter(base_start == min(base_start) | base_end == max(base_end)) %>%
    ungroup()

cna_probes = events %>% filter(event_type == args$cna) %>%
    makeGRangesFromDataFrame(start.field="base_start", end.field="base_end")
probes_overlap = countOverlaps(probes, cna_probes)
n_probes = sum(probes_overlap) / length(probes)
probes$z = scale(probes_overlap)[,1] # max 2.38, sd already << binom approx

top = probes %>%
    seq$intersect(exons) %>%
    group_by_at(vars(-z)) %>%
    summarize(z_aod = mean(z)) %>%
    arrange(-z_aod)

pdf(args$plotfile, 12, 10)
#print(p)
dev.off()
