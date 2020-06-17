library(dplyr)
library(plyranges)
library(ggplot2)
io = import('io')
seq = import('seq')
tcga = import('data/tcga')

plot_compare = function(df, x, y) {
    combined = df %>%
        mutate(cohort = tcga$barcode2study(Sample))

    stats = combined %>%
        na.omit() %>%
        filter(abs({{ x }}) > 0.25 | abs({{ y }}) > 0.25) %>%
        mutate(dist_from_diagonal = abs({{ x }} - {{ y }})) %>%
        group_by(cohort) %>%
        summarize(x = min({{ x }}),
                  y = max({{ y }}),
                  frac_within = sum(dist_from_diagonal < 0.25) / length(dist_from_diagonal))

    ggplot(combined, aes(x={{ x }}, y={{ y }})) +
        stat_bin2d(aes(fill=after_stat(count)), binwidth=c(0.1, 0.1)) +
        scale_fill_distiller(palette="Spectral", trans="log10") +
        geom_abline(slope=1, intercept=c(-0.25,0.25), size=0.5, color="red", linetype="dashed") +
        annotate(geom="rect", xmin=0.75, xmax=1.25, ymin=0.75, ymax=1.25,
                 color="black", linetype="dotted", fill="#ffffff00") +
        geom_text(data=stats, color="red", vjust="inward", hjust="inward",
                  aes(x=x, y=y, label=sprintf("%.0f%% aneup within", 100*frac_within))) +
        facet_wrap(~ cohort, scales="free")
}

cohorts = c("BRCA", "CESC", "ESCA", "LGG", "LUAD", "LUSC", "PRAD", "STAD")

fnames = c("./CESC.WES.Tangent_CBS_TumorsOnly.common_samples_for_IGV.hg19.catted.seg.txt",
           "./LGG_LUSC_PRAD_STAD_252exome_woCNV_hg19.catted.seg.txt")

genes_hg19 = seq$coords$gene(idtype="ensembl_gene_id", assembly="GRCh37", granges=TRUE)
genes_hg38 = seq$coords$gene(idtype="ensembl_gene_id", assembly="GRCh38", granges=TRUE)

exome_copies = lapply(fnames, readr::read_tsv) %>%
    bind_rows() %>%
    mutate(Sample = gsub(".", "-", Sample, fixed=TRUE)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
exome_genes = join_overlap_intersect(genes_hg19, exome_copies) %>%
    as.data.frame() %>%
    as_tibble() %>%
    transmute(Sample = substr(sub("^[^-]+", "TCGA", Sample), 1, 12),
              ensembl_gene_id = ensembl_gene_id,
              exome = log(2^Segment_Mean + 1))

exome2_copies = readr::read_csv("tcga_total_cn_profiles_june1.seg") %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
exome2_genes = join_overlap_intersect(genes_hg19, exome2_copies) %>%
    as.data.frame() %>%
    as_tibble() %>%
    group_by(tcga_barcode) %>%
        mutate(Total_CN = Total_CN / mean(Total_CN, na.rm=TRUE)) %>%
    ungroup() %>%
    transmute(Sample = substr(sub("^[^-]+", "TCGA", tcga_barcode), 1, 12),
              ensembl_gene_id = ensembl_gene_id,
              exome2 = log2(Total_CN + 1))

snp_genes = lapply(cohorts, tcga$cna_genes) %>%
    narray::stack(along=2) %>%
    reshape2::melt() %>%
    as_tibble() %>%
    filter(substr(Var2, 14, 16) == "01A") %>%
    transmute(Sample = substr(gsub(".", "-", Var2, fixed=TRUE), 1, 12),
              ensembl_gene_id = Var1,
              snp = log2(value/2 + 1))

wgs_segments = setdiff(cohorts, c("LUAD", "LUSC")) %>%
    paste0("/data/p282396/data/tcga/TCGAbiolinks-downloader/wgs_segments/TCGA-", ., ".RData") %>%
    io$load() %>%
    bind_rows() %>%
    tcga$filter(along="Sample", cancer=TRUE, primary=TRUE, vial="A") %>%
    mutate(Chromosome = sub("^chr", "", Chromosome)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
wgs_genes = join_overlap_intersect(genes_hg19, wgs_segments) %>%
    as.data.frame() %>%
    as_tibble() %>%
    transmute(Sample = substr(sub("^[^-]+", "TCGA", Sample), 1, 12),
              ensembl_gene_id = ensembl_gene_id,
              wgs = log2(2^Segment_Mean + 1))

ascat_segments = tcga$cna_segments_ascat(granges=TRUE) %>%
    filter(tcga$barcode2study(Sample) %in% cohorts)
ascat_genes = join_overlap_intersect(genes_hg38, ascat_segments) %>%
    as.data.frame() %>%
    as_tibble() %>%
    group_by(Sample) %>%
        mutate(Copy_Number = Copy_Number / mean(Copy_Number, na.rm=TRUE)) %>%
    ungroup() %>%
    transmute(Sample = substr(Sample, 1, 12),
              ensembl_gene_id = ensembl_gene_id,
              ascat = log2(Copy_Number + 1))

pdf("compare.pdf", 16, 10)
print(plot_compare(inner_join(snp_genes, exome_genes), snp, exome))
print(plot_compare(inner_join(snp_genes, exome2_genes), snp, exome2))
print(plot_compare(inner_join(snp_genes, ascat_genes), snp, ascat))
print(plot_compare(inner_join(snp_genes, wgs_genes), snp, wgs))
print(plot_compare(inner_join(wgs_genes, exome_genes), wgs, exome))
print(plot_compare(inner_join(wgs_genes, exome2_genes), wgs, exome2))
print(plot_compare(inner_join(wgs_genes, ascat_genes), wgs, ascat))
dev.off()
