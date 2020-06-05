library(dplyr)
library(plyranges)
library(ggplot2)
io = import('io')
seq = import('seq')
tcga = import('data/tcga')

cohorts = c("BRCA", "CESC", "ESCA", "LGG", "LUAD", "LUSC", "PRAD", "STAD")

fnames = c("./CESC.WES.Tangent_CBS_TumorsOnly.common_samples_for_IGV.hg19.catted.seg.txt",
           "./LGG_LUSC_PRAD_STAD_252exome_woCNV_hg19.catted.seg.txt")

genes_hg19 = seq$coords$gene(idtype="ensembl_gene_id", assembly="GRCh37", granges=TRUE)

exome_copies = lapply(fnames, readr::read_tsv) %>%
    bind_rows() %>%
    mutate(Sample = gsub(".", "-", Sample, fixed=TRUE)) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
exome_genes = join_overlap_intersect(genes_hg19, exome_copies) %>%
    as.data.frame() %>%
    as_tibble() %>%
    transmute(Sample = substr(sub("^[^-]+", "TCGA", Sample), 1, 12),
              ensembl_gene_id = ensembl_gene_id,
              exome = Segment_Mean)

exome2_copies = readr::read_csv("tcga_total_cn_profiles_june1.seg") %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)
exome2_genes = join_overlap_intersect(genes_hg19, exome2_copies) %>%
    as.data.frame() %>%
    as_tibble() %>%
    transmute(Sample = substr(sub("^[^-]+", "TCGA", tcga_barcode), 1, 12),
              ensembl_gene_id = ensembl_gene_id,
              exome2 = log2(Total_CN) - 1)

snp_genes = lapply(cohorts, tcga$cna_genes) %>%
    narray::stack(along=2) %>%
    reshape2::melt() %>%
    as_tibble() %>%
    filter(substr(Var2, 14, 16) == "01A") %>%
    transmute(Sample = substr(gsub(".", "-", Var2, fixed=TRUE), 1, 12),
              ensembl_gene_id = Var1,
              snp = log2(value) - 1)

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
              wgs = Segment_Mean)

#stats = combined %>%
#    na.omit() %>%
#    mutate(dist_from_diagonal = abs(snp - exome)) %>%
#    group_by(cohort) %>%
#    summarize(frac_within = sum(dist_from_diagonal < 0.5) / length(dist_from_diagonal))

pdf("compare.pdf", 16, 10)
#plot(both$snp, both$exome, alpha=0.1)
inner_join(snp_genes, exome_genes) %>%
    mutate(cohort = tcga$barcode2study(Sample)) %>%
    ggplot(aes(x=snp, y=exome)) +
        geom_abline(slope=1, intercept=c(-0.5,0.5), size=0.5, color="red", linetype="dashed") +
        stat_bin2d(aes(fill=after_stat(log(count+1))), binwidth=c(0.1, 0.1)) +
        scale_fill_distiller(palette="Spectral") +
        facet_wrap(~ cohort, scales="free")
#gridExtra::grid.table(stats)

inner_join(snp_genes, exome2_genes) %>%
    mutate(cohort = tcga$barcode2study(Sample)) %>%
    ggplot(aes(x=snp, y=exome2)) +
        geom_abline(slope=1, intercept=c(-0.5,0.5), size=0.5, color="red", linetype="dashed") +
        stat_bin2d(aes(fill=after_stat(log(count+1))), binwidth=c(0.1, 0.1)) +
        scale_fill_distiller(palette="Spectral") +
        facet_wrap(~ cohort, scales="free")

inner_join(snp_genes, wgs_genes) %>%
    mutate(cohort = tcga$barcode2study(Sample)) %>%
    ggplot(aes(x=snp, y=wgs)) +
        geom_abline(slope=1, intercept=c(-0.5,0.5), size=0.5, color="red", linetype="dashed") +
        stat_bin2d(aes(fill=after_stat(log(count+1))), binwidth=c(0.1, 0.1)) +
        scale_fill_distiller(palette="Spectral") +
        facet_wrap(~ cohort, scales="free")

dev.off()
