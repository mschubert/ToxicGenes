library(dplyr)
library(ggplot2)
sys = import('sys')

calc_segs = function(chr) {
    chr %>%
        group_by(bin = round(seq(1, n()/5, length.out=n()))) %>%
        filter(RATIO_CORRECTED > -5) %>%
        summarize(start = min(CHR_POSITION),
                  end = max(CHR_POSITION),
                  mean = mean(RATIO_CORRECTED, na.rm=TRUE, trim=0.2)) %>%
        mutate(seg = ecp::e.divisive(as.matrix(mean), R=20, sig.lvl=0.05)$cluster) %>%
        group_by(seg) %>%
        summarize(start = min(start),
                  end = max(end),
                  mean = mean(mean, na.rm=TRUE))
}

plot_cnv = function(cnv) {
    pt = unique(cnv$PtID)
    dset = inner_join(cnv, meta)
    segs = segs %>% filter(PtID == pt)
    ggplot(dset, aes(x=CHR_POSITION, y=RATIO_CORRECTED)) +
        geom_point(size=0.2, alpha=0.4) +
        geom_segment(data=segs, aes(x=start, xend=end, y=mean, yend=mean), color="firebrick") +
        facet_grid(. ~ chr, scales="free_x", space="free_x") +
        ggtitle(pt)
}

args = sys$cmd$parse(
    opt('i', 'indir', 'txt dir', 'N4plus_CNVseq_Michael_20230414'),
    opt('m', 'meta', 'xlsx', 'N4plus_data_Michael_20230414.xlsx'),
    opt('o', 'outfile', 'rds', 'plot_cnvs.rds'),
    opt('p', 'plotfile', 'pdf', 'plot_cnvs.pdf')
)

meta = readxl::read_xlsx(args$meta) %>%
    mutate_all(function(x) ifelse(x == 9999, NA, x)) %>%
    mutate(PtID = sprintf("%03i", PtID))
fnames = gtools::mixedsort(list.files(args$indir, full.names=TRUE))
cnvs = setNames(fnames, sub("^([0-9]+)-.*", "\\1", basename(fnames))) %>%
    lapply(readr::read_tsv) %>%
    bind_rows(.id="PtID") %>%
    filter(CHROMOSOME %in% paste0("chr", c(1:22,'X'))) %>%
    mutate(chr = factor(CHROMOSOME, levels=gtools::mixedsort(unique(CHROMOSOME))))

segs = cnvs %>% group_by(PtID, chr) %>%
    tidyr::nest() %>%
    ungroup() %>%
    mutate(segs = clustermq::Q(calc_segs, chr=data, n_jobs=50, pkgs="dplyr")) %>%
    select(-data) %>%
    tidyr::unnest(cols=segs)

saveRDS(list(cnvs=cnvs, segs=segs), file=args$outfile)

plots = split(cnvs, cnvs$PtID) %>% lapply(plot_cnv)
pdf(args$plotfile, 12, 5)
for (p in plots)
    print(p)
dev.off()
