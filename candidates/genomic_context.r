library(ggplot2)
library(patchwork)
library(dplyr)
library(plyranges)
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga')

plot_gene_track = function(gene="CDKN1A", et=0.15, len=1e8, bins=1000, hl=len/100) {
    pos = as.data.frame(genes_hg38) %>%
        filter(external_gene_name == gene)
    center_pos = round(mean(with(pos, c(start, end))))

    dset = data_frame(seqnames = pos$seqnames,
                      start=seq(1, len, by=len/bins) + max(0, center_pos - len/2)) %>%
        mutate(end = start + len/bins - 1) %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        join_overlap_inner_within(copies %>% select(Sample, ploidy)) %>%
        group_by(start, end) %>%
        summarize(eup = n_distinct(Sample),
                  gain = n_distinct(Sample[ploidy>2+et & ploidy<3-et]),
                  amp = n_distinct(Sample[ploidy>3-et]),
                  loss = -n_distinct(Sample[ploidy<2-et & ploidy>1+et]),
                  del = -n_distinct(Sample[ploidy<1+et])) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(gain = 100*gain/eup, amp=100*amp/eup, loss=100*loss/eup, del=100*del/eup) %>%
        select(-eup) %>%
        tidyr::gather("type", "num", -start, -end)

    genes = data.frame(seqnames = pos$seqnames,
                       start = max(0, center_pos-len/2)) %>%
        mutate(end = start + len) %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        join_overlap_intersect(genes_hg38 %>% select(external_gene_name)) %>%
        as.data.frame() %>%
        as_tibble()

    cols = c(gain="tomato", amp="firebrick", loss="lightblue", del="blue")
    p1 = ggplot(dset, aes(x=start)) +
        geom_vline(xintercept=c(center_pos-hl/2, center_pos+hl/2), linetype="dotted") +
        geom_ribbon(aes(ymax=num, group=type, fill=type), ymin=0, alpha=0.2) +
        geom_line(aes(y=num, color=type)) +
        scale_fill_manual(values=cols, guide=FALSE) +
        scale_color_manual(values=cols, guide=FALSE) +
        theme_classic() +
        geom_segment(data=genes, aes(xend=end), y=0, yend=0, size=2, alpha=0.7) +
        ggtitle(sprintf("%s chr%i:%.1f-%.1f Mb", gene, pos$seqnames, min(dset$start)/1e6, max(dset$end)/1e6)) +
        theme(axis.title.x = element_blank())

    p2 = p1 +
        xlim(c(center_pos-hl/2, center_pos+hl/2)) +
        ggrepel::geom_text_repel(data=genes, size=2, max.iter=1e5, segment.alpha=0.3,
            aes(x=0.5*(start+end), y=0, label=external_gene_name)) +
        ggtitle(sprintf("%.1f-%.1f Mb", (center_pos-hl/2)/1e6, (center_pos+hl/2)/1e6)) +
        theme(axis.title.x = element_blank())

    p1 / p2
}

args = sys$cmd$parse(
    opt('c', 'cohort', 'chr', 'pan'),
#    opt('', '', '', ''),
    opt('p', 'plotfile', 'pdf', 'gctx_pan.pdf')
)

if (args$cohort == "pan")
    args$cohort = c("BRCA", "LUAD", "LUSC", "OV", "PRAD", "SKCM")

genes_hg38 = seq$coords$gene(idtype="hgnc_symbol", assembly="GRCh38", granges=TRUE)
purity = tcga$purity() %>%
    select(Sample, purity=estimate) %>%
    na.omit()
copies = lapply(args$cohort, tcga$cna_segments) %>%
    bind_rows() %>%
    inner_join(purity) %>%
    mutate(ploidy = (ploidy-2)/purity + 2) %>% #todo: if purity/SegMean incorrect that can be <0
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
#racs =

pdf(args$plotfile, 10, 4)
for (g in c("EGFR", "ERBB2", "MYC", "CDKN1A", "RBM14", "RBM12", "HES1",
            "SNRPA", "NSF", "BANP", "H3F3C", "IRF2", "ZBTB14"))
    print(plot_gene_track(g))
dev.off()
