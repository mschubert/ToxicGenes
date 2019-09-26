library(dplyr)
library(ggplot2)

events = readr::read_tsv("./TCGA.pancan12.130815.zigg_events.160923.txt") 

real_ev = events %>%
    filter(amplitude >= 0.1,
           ! (arm_length > 0.95 & arm_length < 1.05),
           ! (arm_length > 1.95))

plot_chr = function(chrom) {
    message(chrom)
    chr = real_ev %>%
        filter(chr == chrom)

    pos = sample(1e7:(max(chr$base_end)-1e7), 1e4, replace=TRUE)

    pos2left = function(p) { #todo: filter out centromere spanning?
        cur = chr %>%
            filter(base_start < p & base_end > p)
        re = data.frame(left=p - cur$base_start, right = cur$base_end - p)
    }
    resmp = lapply(pos, pos2left) %>%
        dplyr::bind_rows()
    p = resmp %>%
        tidyr::gather("side", "end") %>%
        ggplot(aes(x=end)) +
            geom_histogram(bins=100) +
    #        stat_smooth(method=nls, formula=..count.. ~ 1/x) +
            facet_wrap(~ side) +
            ggtitle(sprintf("chr%s random positions", chrom))
}

pdf("ends.pdf", 12, 8)
for (chr in c(1:22,'X'))
    print(plot_chr(chr))
dev.off()
