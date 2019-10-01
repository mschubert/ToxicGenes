library(dplyr)
library(ggplot2)

events = readr::read_tsv("./TCGA.pancan12.130815.zigg_events.160923.txt") 

p = ggplot(events, aes(x=base_end - base_start, y=arm_length)) +
    geom_point() +
    facet_wrap(~ chr)

png("cor_arm_length.png", 2000, 1800)
print(p)
dev.off()
