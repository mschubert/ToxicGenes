library(ggplot2)
library(patchwork)
library(dplyr)
library(plyranges)
sys = import('sys')
seq = import('seq')
gdsc = import('data/gdsc')
tcga = import('data/tcga')

rmean_noNA = function(x, size=3) {
    nona = ! is.na(x)
    x[nona] = zoo::rollapply(x[nona], size, mean, partial=TRUE, fill=0)
    x
}

# needs objs: genes_hg38, racs, copies
plot_gene_track = function(gene="CDKN1A", et=0.15, len=1e8, bins=1000, hl=len/100) {
    message(gene)
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
        as.data.frame() %>% as_tibble() %>%
        select(-eup) %>%
        tidyr::gather("type", "num", -start, -end)

    genes = data.frame(seqnames = pos$seqnames,
                       start = max(0, center_pos-len/2)) %>%
        mutate(end = start + len) %>%
        GenomicRanges::makeGRangesFromDataFrame() %>%
        join_overlap_intersect(genes_hg38 %>% select(external_gene_name, driver, comp_tcga_rmean, comp_ccle_rmean)) %>%
        as.data.frame() %>%
        as_tibble() %>%
        mutate(midpoint = (start+end)/2,
               external_gene_name = ifelse(end < center_pos-hl/2 | start > center_pos+hl/2,
                                           NA, external_gene_name))
    comp = genes %>%
        select(seqnames, midpoint, comp_tcga_rmean, comp_ccle_rmean) %>%
        tidyr::gather("type", "rmean", -seqnames, -midpoint) %>%
        na.omit()

    arws = c(Amplification="▲", Deletion="▼")
    ys = c(Amplification=max(dset$num)/2, Deletion=min(dset$num)/2)
    cracs = tibble(seqnames = pos$seqnames,
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
    cna_range = diff(range(dset$num))
    comp_mult = cna_range / 2
    p0 = ggplot(dset, aes(x=start)) +
#        geom_rect(data=fras, aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, fill="grey35", alpha=0.1) +
#        geom_rect(data=cracs, aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf, fill="gold", alpha=0.1) +
        geom_vline(xintercept=c(center_pos-hl/2, center_pos+hl/2), linetype="dotted") +
        geom_ribbon(aes(ymax=num, group=type, fill=type), ymin=0, alpha=0.2) +
        geom_line(aes(y=num, color=type)) +
        geom_segment(data=genes, aes(xend=end), y=0, yend=0, size=2, alpha=0.5, color="black") +
        geom_segment(data=genes %>% filter(external_gene_name == gene), aes(xend=end),
                     y=0, yend=0, size=2, alpha=1, color="deeppink") +
        geom_density(data=genes, aes(y=..scaled..*(cna_range/5) + (cna_range/50)),
                     alpha=0.1, fill="#ffffff00", bw=1.5e5, size=0.4) +
        geom_line(data=comp %>% filter(type=="comp_tcga_rmean"), aes(x=midpoint, y=rmean*comp_mult),
                  color="#7cae00", size=0.2) +
        geom_line(data=comp %>% filter(type=="comp_ccle_rmean"), aes(x=midpoint, y=rmean*comp_mult),
                  color="#c77cff", size=0.2) +
        scale_fill_manual(values=cols, guide=FALSE) +
        scale_color_manual(values=cols, guide=FALSE) +
        scale_y_continuous(name="# samples", sec.axis = sec_axis(~./comp_mult, name="compensation")) +
        theme_classic() +
        theme(axis.title.x = element_blank())

    p1 = p0 +
#        ggrepel::geom_text_repel(data=labs, aes(x=x, y=y, label=label), size=2, hjust=0.5) +
        ggrepel::geom_label_repel(data=genes, size=2, max.iter=1e5, segment.alpha=0.5,
            aes(x=0.5*(start+end), y=0, label=driver), min.segment.length=0,
            label.size=NA, fill="#ffffff80", segment.color="white") +
        coord_cartesian(expand=FALSE) +
        ggtitle(sprintf("%s chr%i:%.0f-%.0f Mb", gene, pos$seqnames, min(dset$start)/1e6, max(dset$end)/1e6))

    p2 = p0 +
        coord_cartesian(xlim=c(center_pos-hl/2, center_pos+hl/2)) +
        ggrepel::geom_label_repel(data=genes, size=2, max.iter=1e5, segment.alpha=0.5,
            aes(x=0.5*(start+end), y=0, label=external_gene_name),
            min.segment.length=0, label.size=NA, fill="#ffffff80", segment.color="white") +
        ggtitle(sprintf("%.1f-%.1f Mb", (center_pos-hl/2)/1e6, (center_pos+hl/2)/1e6))

    p1 / p2
}

args = sys$cmd$parse(
    opt('c', 'cohort', 'chr', 'pan'),
    opt('m', 'comp_tcga', 'tcga xlsx', '../tcga/pan/stan-nb_puradj.xlsx'),
    opt('n', 'comp_ccle', 'ccle xlsx', '../ccle/pan/stan-nb.xlsx'),
    opt('p', 'plotfile', 'pdf', 'gctx_pan.pdf')
)

if (args$cohort == "pan") {
    args$cohort = c("BRCA", "LUAD", "LUSC", "OV", "PRAD", "SKCM")
    racs_cohort = "PANCAN"
} else if (args$cohort == "COADREAD") {
    args$cohort = c("COAD", "READ")
    racs_cohort = "COAD/READ"
} else if (args$cohort == "NSCLC") {
    racs_cohort = args$cohort = c("LUAD", "LUSC")
} else {
    racs_cohort = args$cohort
}

tcga_comp = readxl::read_xlsx(args$comp_tcga) %>% #, "amp") %>%
    transmute(external_gene_name = gene, comp_tcga = pmax(pmin(estimate, 3), -3))
ccle_comp = readxl::read_xlsx(args$comp_ccle) %>%
    transmute(external_gene_name = gene, comp_ccle = pmax(pmin(estimate, 3), -3))
fra = readr::read_tsv("../data/fragile_sites/fra.txt") %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
seqlevelsStyle(fra) = "Ensembl"
drivers = gdsc$drivers(args$cohort) %>% pull(HGNC) %>% unique()
genes_hg38 = seq$coords$gene(idtype="hgnc_symbol", assembly="GRCh38") %>%
    filter(gene_biotype %in% c("protein_coding", "miRNA")) %>%
    mutate(driver = ifelse(external_gene_name %in% drivers, external_gene_name, NA)) %>%
    left_join(tcga_comp) %>%
    left_join(ccle_comp) %>%
    mutate(strand = setNames(c("+", "-"), c(1, -1))[as.character(strand)]) %>%
    GenomicRanges::makeGRangesFromDataFrame(start.field="start_position",
            end.field="end_position", keep.extra.columns=TRUE) %>%
    arrange(seqnames, (start+end)/2) %>%
    group_by(seqnames) %>%
        mutate(comp_tcga_rmean = rmean_noNA(comp_tcga),
               comp_ccle_rmean = rmean_noNA(comp_ccle)) %>%
    ungroup()
purity = tcga$purity() %>%
    select(Sample, purity=consensus) %>%
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
            "SNRPA", "NSF", "BANP", "IRF2", "ZBTB14"))
    print(plot_gene_track(g))
dev.off()
