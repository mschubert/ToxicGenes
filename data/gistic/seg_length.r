library(dplyr)
library(ggplot2)

events = readr::read_tsv("./TCGA.pancan12.130815.zigg_events.160923.txt") 

real_ev = events %>%
    filter(amplitude >= 0.1,
           arm_length < 1.95,
           arm_length < 0.95 | arm_length > 1.05)

plot_chr = function(chrom) {
    message(chrom)
    chr = real_ev %>%
        mutate(len = base_end - base_start) %>%
        filter(chr == chrom,
               len <= 1e8)

    p = chr %>%
        ggplot(aes(x=len)) +
            geom_histogram(bins=100) +
            scale_y_log10() +
    #        stat_smooth(method=nls, formula=..count.. ~ 1/x) +
            ggtitle(sprintf("chr%s segment lengths", chrom))
}

pdf("seg_length.pdf", 10, 8)
for (chr in c(1:22,'X'))
    print(plot_chr(chr))
dev.off()
