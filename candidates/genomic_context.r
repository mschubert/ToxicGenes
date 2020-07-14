library(ggplot2)
library(patchwork)
library(dplyr)
library(plyranges)
sys = import('sys')
seq = import('seq')
gdsc = import('data/gdsc')
tcga = import('data/tcga')

rmean_noNA = function(x) {
    nona = ! is.na(x)
    x[nona] = zoo::rollapply(x[nona], 5, mean, partial=TRUE, fill=0)
    x
}

# needs objs: genes_hg38, racs, copies
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
        join_overlap_intersect(genes_hg38 %>% select(external_gene_name, driver, comp_pct, comp_pct_rmean)) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(external_gene_name = ifelse(end < center_pos-hl/2 | start > center_pos+hl/2,
                                           NA, external_gene_name))

    comp = genes %>%
        filter(!is.na(comp_pct_rmean)) %>%
        transmute(midpoint = (start+end)/2, comp_pct, comp_pct_rmean)

    arws = c(Amplification="▲", Deletion="▼")
    ys = c(Amplification=max(dset$num)/2, Deletion=min(dset$num)/2)
    cracs = data_frame(seqnames = pos$seqnames,
                       start=max(0, center_pos - len/2)) %>%
        mutate(end = start + len) %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        join_overlap_inner_within(racs %>% select(racs, cna), .) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(label=paste(racs, arws[cna]), x=(start+end)/2, y=ys[cna])

    fras = data_frame(seqnames = pos$seqnames,
                      start=max(0, center_pos - len/2)) %>%
        mutate(end = start + len) %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        join_overlap_inner_within(fra, .) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(label=fra, x=(start+end)/2, y = max(dset$num)/2)

    labs = bind_rows(cracs %>% select(x, y, label),
                     fras %>% select(x, y, label)) %>% na.omit()

    cols = c(gain="tomato", amp="firebrick", loss="lightblue", del="blue")
    p0 = ggplot(dset, aes(x=start)) +
        geom_rect(data=fras, aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, fill="grey35", alpha=0.1) +
        geom_rect(data=cracs, aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, fill="gold", alpha=0.1) +
        geom_vline(xintercept=c(center_pos-hl/2, center_pos+hl/2), linetype="dotted") +
        geom_ribbon(aes(ymax=num, group=type, fill=type), ymin=0, alpha=0.2) +
        geom_line(aes(y=num, color=type)) +
        scale_fill_manual(values=cols, guide=FALSE) +
        scale_color_manual(values=cols, guide=FALSE) +
        theme_classic() +
        geom_segment(data=genes, aes(xend=end), y=0, yend=0, size=2, alpha=0.5) +
        theme(axis.title.x = element_blank())

    p1 = p0 +
        geom_line(data=comp, aes(x=midpoint, y=comp_pct_rmean)) +
        ggrepel::geom_text_repel(data=labs, aes(x=x, y=y, label=label), size=2, hjust=0.5) +
        ggrepel::geom_label_repel(data=genes, size=2, max.iter=1e5, segment.alpha=0.5,
            aes(x=0.5*(start+end), y=0, label=driver), min.segment.length=0,
            label.size=NA, fill="#ffffff80", segment.color="white") +
        ggtitle(sprintf("%s chr%i:%.0f-%.0f Mb", gene, pos$seqnames, min(dset$start)/1e6, max(dset$end)/1e6))

    p2 = p0 +
        geom_line(data=comp, aes(x=midpoint, y=comp_pct)) +
        coord_cartesian(xlim=c(center_pos-hl/2, center_pos+hl/2)) +
        ggrepel::geom_label_repel(data=genes, size=2, max.iter=1e5, segment.alpha=0.5,
            aes(x=0.5*(start+end), y=0, label=external_gene_name),
            min.segment.length=0, label.size=NA, fill="#ffffff80", segment.color="white") +
        ggtitle(sprintf("%.1f-%.1f Mb", (center_pos-hl/2)/1e6, (center_pos+hl/2)/1e6))

    p1 / p2
}

args = sys$cmd$parse(
    opt('c', 'cohort', 'chr', 'pan'),
    opt('m', 'compensation', 'tcga xlsx', '../tcga/pan/rlm3_naive.xlsx'),
    opt('p', 'plotfile', 'pdf', 'gctx_pan.pdf')
)

if (args$cohort == "pan") {
    args$cohort = c("BRCA", "LUAD", "LUSC", "OV", "PRAD", "SKCM")
    racs_cohort = "PANCAN"
} else {
    racs_cohort = args$cohort
}

tcga_comp = readxl::read_xlsx(args$compensation, "amp") %>%
    transmute(external_gene_name = name,
              comp_pct = statistic)
fra = readr::read_tsv("../data/fragile_sites/fra.txt") %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
seqlevelsStyle(fra) = "Ensembl"
drivers = gdsc$drivers(args$cohort) %>% pull(HGNC) %>% unique()
genes_hg38 = seq$coords$gene(idtype="hgnc_symbol", assembly="GRCh38") %>%
    filter(gene_biotype %in% c("protein_coding", "miRNA")) %>%
    mutate(driver = ifelse(external_gene_name %in% drivers, external_gene_name, NA)) %>%
    left_join(tcga_comp) %>%
    mutate(strand = setNames(c("+", "-"), c(1, -1))[as.character(strand)]) %>%
    GenomicRanges::makeGRangesFromDataFrame(start.field="start_position",
            end.field="end_position", keep.extra.columns=TRUE) %>%
    arrange(seqnames, (start+end)/2) %>%
    group_by(seqnames) %>%
        mutate(comp_pct_rmean = rmean_noNA(comp_pct)) %>%
    ungroup()
purity = tcga$purity() %>%
    select(Sample, purity=estimate) %>%
    na.omit()
copies = lapply(args$cohort, tcga$cna_segments) %>%
    bind_rows() %>%
    inner_join(purity) %>%
    mutate(ploidy = (ploidy-2)/purity + 2) %>% #todo: if purity/SegMean incorrect that can be <0
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
hg19_to_hg38 = rtracklayer::import.chain(system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain"))
racs = readxl::read_xlsx("../data/racs/TableS2D.xlsx", skip=19) %>%
    setNames(make.names(colnames(.))) %>%
    select(racs=Identifier, cohort=Cancer.Type, chr, start, stop,
           cna=Recurrent..Amplification...Deletion, Contained.genes) %>%
    filter(cohort == racs_cohort) %>%
    select(-cohort, -Contained.genes) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
seqlevelsStyle(racs) = "UCSC"
racs2 = rtracklayer::liftOver(racs, hg19_to_hg38) %>%
    unlist() %>%
    group_by(seqnames, racs, cna) %>%
    summarize(start=min(start), end=max(end)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
seqlevelsStyle(racs2) = "Ensembl"
racs = racs2

cairo_pdf(args$plotfile, 10, 4, onefile=TRUE)
for (g in c("EGFR", "ERBB2", "MYC", "CDKN1A", "RBM14", "RBM12", "HES1",
            "SNRPA", "NSF", "BANP", "H3F3C", "IRF2", "ZBTB14"))
    print(plot_gene_track(g))
dev.off()
