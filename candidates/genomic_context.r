library(ggplot2)
library(dplyr)
library(plyranges)
seq = import('seq')
tcga = import('data/tcga')

et = 0.15
len = 1e8
bins = 1000
hl = len/100

genes_hg38 = seq$coords$gene(idtype="hgnc_symbol", assembly="GRCh38", granges=TRUE)
copies = tcga$cna_segments("BRCA", granges=TRUE) # test on BRCA
#racs =

center = "CDKN1A"
pos = as.data.frame(genes_hg38) %>%
    filter(external_gene_name == center)
center_pos = round(mean(with(pos, c(start, end))))

dset = data_frame(seqnames = pos$seqnames,
                  start=seq(1, len, by=len/bins) + max(0, center_pos - len/2)) %>%
    mutate(end = start + len/bins - 1) %>%
    GenomicRanges::makeGRangesFromDataFrame() %>%
    join_overlap_inner_within(copies %>% select(Sample, ploidy)) %>%
    group_by(start, end) %>%
    summarize(eup = n_distinct(Sample[abs(ploidy-2)<et]),
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
p = ggplot(dset, aes(x=start)) +
    geom_rect(aes(xmin=center_pos-hl/2, xmax=center_pos+hl/2), ymin=-Inf, ymax=Inf, fill="lightgrey") +
    geom_ribbon(aes(ymax=num, group=type, fill=type), ymin=0, alpha=0.2) +
    geom_line(aes(y=num, color=type)) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) +
    theme_classic() +
    geom_segment(data=genes, aes(xend=end), y=0, yend=0, size=2, alpha=0.7) +
#    ggrepel::geom_text_repel(data=genes, size=2, maxiter=1e5, segment.alpha=0.3,
#        aes(x=0.5*(start+end), y=0, label=external_gene_name)) +
    ggtitle(sprintf("%s chr%i:%.1f-%.1f Mb", center, pos$seqnames, min(dset$start)/1e6, max(dset$end)/1e6))

pdf("genomic_context.pdf", 10, 4)
print(p)
dev.off()
