library(dplyr)
library(plyranges)
library(ggplot2)
seq = import('seq')
tcga = import('data/tcga')

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

snp_genes = c("CESC", "LGG", "LUSC", "PRAD", "STAD") %>%
    lapply(tcga$cna_genes) %>%
    narray::stack(along=2) %>%
    reshape2::melt() %>%
    as_tibble() %>%
    filter(substr(Var2, 14, 16) == "01A") %>%
    transmute(Sample = substr(gsub(".", "-", Var2, fixed=TRUE), 1, 12),
              ensembl_gene_id = Var1,
              snp = log2(value) - 1)

both = inner_join(snp_genes, exome_genes) %>%
    mutate(cohort = tcga$barcode2study(Sample)) %>%
    filter(abs(snp) + abs(exome) > 0.5)

stats = both %>%
    na.omit() %>%
    mutate(dist_from_diagonal = abs(snp - exome)) %>%
    group_by(cohort) %>%
    summarize(frac_within = sum(dist_from_diagonal < 0.5) / length(dist_from_diagonal))

pdf("compare.pdf", 16, 10)
#plot(both$snp, both$exome, alpha=0.1)
ggplot(both, aes(x=snp, y=exome)) +
    geom_abline(slope=1, intercept=c(-0.5,0.5), size=0.5, color="red", linetype="dashed") +
    stat_bin2d(aes(fill=after_stat(log(count+1))), binwidth=c(0.1, 0.1)) +
    scale_fill_distiller(palette="Spectral") +
    facet_wrap(~ cohort, scales="free")
gridExtra::grid.table(stats)
dev.off()
