library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'plot_cnvs.rds'),
    opt('m', 'meta', 'xlsx', 'N4plus_data_Michael_20230414.xlsx'),
    opt('p', 'plotfile', 'pdf', 'surv.pdf')
)

meta = readxl::read_xlsx(args$meta) %>%
    mutate_all(function(x) ifelse(x == 9999, NA, x)) %>%
    mutate(PtID = sprintf("%03i", PtID))

dset = readRDS(args$infile)
segs = dset$segs %>%
    filter(chr == "chr11")

genes = tibble(genes = c("CCND1", "RBM14"), pos=c(69647000, 66622000))

p1 = ggplot(segs) +
    geom_hline(yintercept=0, size=1, color="red") +
    geom_segment(aes(x=start, xend=end, y=mean, yend=mean), alpha=0.1) +
    geom_point(data=genes, aes(x=pos), y=0, color="red", size=3, alpha=0.5) +
    ggrepel::geom_label_repel(data=genes, aes(x=pos, label=genes), y=0,
                              color="red", alpha=0.5) +
    coord_cartesian(expand=FALSE)

segs2 = segs %>%
    filter(start < 69647000, end > 66622000)
p2 = ggplot(segs2) +
    geom_hline(yintercept=0, size=1, color="red") +
    geom_vline(data=genes, aes(xintercept=pos), linetype="dashed", color="red", alpha=0.2) +
    geom_segment(aes(x=start, xend=end, y=mean, yend=mean), alpha=0.1) +
    geom_point(data=genes, aes(x=pos), y=0, color="red", size=3, alpha=0.5) +
    ggrepel::geom_label_repel(data=genes, aes(x=pos, label=genes), y=0,
                              color="red", alpha=0.5) +
    coord_cartesian(expand=FALSE)

pdf(args$plotfile)
print(p1)
print(p2)
dev.off()
