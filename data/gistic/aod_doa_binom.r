library(dplyr)
library(ggplot2)
library(GenomicRanges)

events = readr::read_tsv("./TCGA.pancan12.130815.zigg_events.160923.txt")
probes = readr::read_tsv("genome.info.6.0_hg19.na31_minus_frequent_nan_probes_sorted_2.1 2.txt",
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

aod = events %>% filter(event_type == "aod") %>%
    makeGRangesFromDataFrame(start.field="base_start", end.field="base_end")
probes_aod = countOverlaps(probes, aod)
np_aod = sum(probes_aod) / length(probes)
z = scale(probes_aod) # max 2.38, sd already << binom approx

# binom.test(0, n=np_aod, p=1/l..)
# N(np_aod/length(probes), )

doa = events %>% filter(event_type == "doa") %>%
    makeGRangesFromDataFrame(start.field="base_start", end.field="base_end")
probes_doa = countOverlaps(probes, doa)
np_doa = sum(probes_doa) / length(probes)
z = scale(probes_doa) # max 3.09 (p=0.002)



pdf("aod_doa_binom.pdf", 12, 10)
print(p)
dev.off()
