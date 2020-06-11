library(dplyr)
library(plyranges)
library(ggplot2)
io = import('io')
seq = import('seq')
ccle = import('data/ccle')

plot_compare = function(combined, x, y) {
    stats = combined %>%
        na.omit() %>%
        filter(abs({{ x }}-1) > 0.25 | abs({{ y }}-1) > 0.25) %>%
        mutate(dist_from_diagonal = abs({{ x }} - {{ y }})) %>%
#        group_by(primary_disease) %>%
        summarize(x = min({{ x }}),
                  y = max({{ y }}),
                  frac_within = sum(dist_from_diagonal < 0.25) / length(dist_from_diagonal))

    ggplot(combined, aes(x={{ x }}, y={{ y }})) +
        stat_bin2d(aes(fill=after_stat(log(count+1))), binwidth=c(0.1, 0.1)) +
        scale_fill_distiller(palette="Spectral") +
        geom_abline(slope=1, intercept=c(-0.25,0.25), size=0.5, color="red", linetype="dashed") +
        annotate(geom="rect", xmin=0.75, xmax=1.25, ymin=0.75, ymax=1.25,
                 color="black", linetype="dotted", fill="#ffffff00") +
        geom_text(data=stats, color="red", vjust="inward", hjust="inward",
                  aes(x=x, y=y, label=sprintf("%.0f%% aneup within", 100*frac_within))) +
#        facet_wrap(~ primary_disease, scales="free")
        coord_fixed()
}

genes = seq$coords$gene(idtype="hgnc_symbol", assembly="GRCh38", granges=TRUE)

info = readr::read_csv("sample_info.csv")
exome_genes = readr::read_csv("./CCLE_gene_cn.csv") %>%
    tidyr::gather("gene", "exome", -X1) %>%
    dplyr::rename(DepMap_ID=X1) %>%
    mutate(gene = sub(" .[0-9]+.$", "", gene))

#todo: check if this is GRCh38
absolute_segments = readxl::read_xlsx("CCLE_ABSOLUTE_combined_20181227.xlsx", na="NA") %>%
    dplyr::rename(CCLE_Name=sample) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=TRUE)
absolute_genes = join_overlap_intersect(genes, absolute_segments) %>%
    as.data.frame() %>% as_tibble() %>%
    transmute(DepMap_ID = depMapID,
              gene = external_gene_name,
              absolute = Modal_Total_CN + Subclonal_HSCN_a1 + Subclonal_HSCN_a2) %>%
    mutate(absolute = log2(absolute/2+1)) # same as exome

snp_genes = readRDS("../ccle/dset.rds")$copies %>%
    reshape2::melt() %>%
    dplyr::rename(gene=Var1, CCLE_Name=Var2, snp=value) %>%
    as_tibble() %>%
    mutate(snp = log2(snp/2+1)) # same as exome

combined = info %>%
    select(DepMap_ID, CCLE_Name, primary_disease) %>% # lineage?
    tidyr::crossing(genes %>% as.data.frame() %>% select(gene=external_gene_name)) %>%
    left_join(snp_genes) %>%
    left_join(exome_genes) %>%
    left_join(absolute_genes) %>%
    filter(is.na(snp)+ is.na(exome) + is.na(absolute) < 2)

pdf("compare.pdf", 10, 8)
print(plot_compare(combined, snp, exome))
print(plot_compare(combined, snp, absolute))
print(plot_compare(combined, exome, absolute))
dev.off()
