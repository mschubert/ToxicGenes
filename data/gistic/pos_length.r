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
        filter(chr == chrom)

    pos = sample(2*1e7:(max(chr$base_end)-2*1e7), 1e4, replace=TRUE)

    pos2left = function(p) { #todo: filter out centromere spanning?
        cur = chr %>%
            filter(base_start < p & base_end > p)
        re = data.frame(len = cur$base_end - cur$base_start)
    }
    resmp = lapply(pos, pos2left) %>%
        dplyr::bind_rows()
    p = resmp %>%
        filter(len < 1e7) %>%
        ggplot(aes(x=len)) +
            geom_histogram(bins=100) +
    #        stat_smooth(method=nls, formula=..count.. ~ 1/x) +
            ggtitle(sprintf("chr%s random positions lengths", chrom))
}

pdf("pos_length.pdf", 10, 8)
for (chr in c(1:22,'X'))
    print(plot_chr(chr))
dev.off()
