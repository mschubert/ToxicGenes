library(ggplot2)
library(dplyr)
library(plyranges)
seq = import('seq')
tcga = import('data/tcga')

genes_hg38 = seq$coords$gene(idtype="hgnc_symbol", assembly="GRCh38", granges=TRUE)
copies = tcga$cna_segments("BRCA", granges=TRUE) # test on BRCA
#racs =

center = "CDKN1A"
pos = as.data.frame(genes_hg38) %>%
    filter(external_gene_name == center)
center_pos = round(mean(with(pos, c(start, end))))

et = 0.15

len = 1e7
bins = 100
dset = data_frame(seqnames = pos$seqnames,
                  start=seq(1, len, by=len/bins) + center_pos - len/2) %>%
    mutate(end = start + len/bins - 1) %>%
    GenomicRanges::makeGRangesFromDataFrame() %>%
    join_overlap_left_within(copies %>% select(Sample, ploidy)) %>%
    group_by(start, end) %>%
    summarize(eup = sum(abs(ploidy-2)<et),
              gain = sum(ploidy>2+et & ploidy<3-et),
              amp = sum(ploidy>3-et),
              loss = -sum(ploidy<2-et & ploidy>1+et),
              del = -sum(ploidy<1+et)) %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(gain = 100*gain/eup, amp=100*amp/eup, loss=100*loss/eup, del=100*del/eup) %>%
    select(-eup) %>%
    tidyr::gather("type", "num", -start, -end)

genes = data.frame(seqnames=pos$seqnames, start=center_pos-len/2, end=center_pos+len/2) %>%
    GenomicRanges::makeGRangesFromDataFrame() %>%
    join_overlap_intersect(genes_hg38 %>% select(external_gene_name)) %>%
    as.data.frame() %>%
    as_tibble()

cols = c(gain="tomato", amp="firebrick", loss="lightblue", del="blue")
ggplot(dset, aes(x=start)) +
    geom_ribbon(aes(ymax=num, group=type, fill=type), ymin=0, alpha=0.2) +
    geom_line(aes(y=num, color=type)) +
    scale_fill_manual(values=cols) +
    scale_color_manual(values=cols) +
    theme_classic() +
    geom_segment(data=genes, aes(xend=end), y=0, yend=0, size=2, alpha=0.7) +
#    ggrepel::geom_text_repel(data=genes, size=2, maxiter=1e5, segment.alpha=0.3,
#        aes(x=0.5*(start+end), y=0, label=external_gene_name)) +
    ggtitle(sprintf("%s chr%i:%i-%i", center, pos$seqnames, min(dset$start), max(dset$end)))
