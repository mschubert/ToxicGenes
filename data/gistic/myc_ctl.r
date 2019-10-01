library(dplyr)
library(ggplot2)

#TODO:
# filter out telomere-bound events
# comp.: BrISCUT (as unexpected cutoff)
# comp. as DOA enrichment over D? (also: absence of amp/del, but need to adjust that for drivers?)

filter_plot = function(ev, ch) {
    myc = ev %>%
        filter(chr == ch, amplitude >= 0.1,
               base_start <= myc_start & base_end >= myc_end)

    ggplot(myc, aes(x=base_end - base_start)) +
        geom_histogram(bins=100)
}

# hg19 c(start, end)
myc_start = 128747680
myc_end = 128753680

# myc fully contained in CNA
events = readr::read_tsv("./TCGA.pancan12.130815.zigg_events.160923.txt") 

pdf("myc.pdf", 12, 6)

ggplot(events, aes(x=base_end - base_start)) +
    geom_histogram(bins=100) +
    ggtitle("everything")

filter_plot(events, 8) +
    ggtitle(sprintf("covering MYC (chr8:%i-%i)", myc_start, myc_end))

for (i in c(1:7)) {
    p = filter_plot(events, i) +
        ggtitle(sprintf("chr%i:%i-%i", i, myc_start, myc_end))
    print(p)
}

dev.off()
