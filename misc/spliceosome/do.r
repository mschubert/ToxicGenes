library(dplyr)
library(ggplot2)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')

# subset splicosome
splice = readr::read_tsv("spliceosome.txt") %>%
    select(ensembl_gene_id = `Ensembl gene ID`,
           hgnc_symbol = `Approved symbol`,
           chr = Chromosome,
           subunit = `Group name`) %>%
    distinct()

# load breast cancer_copies
purity = tcga$purity() %>%
    filter(!is.na(estimate))
copies = tcga$cna_genes("BRCA")
narray::intersect(purity$Sample, copies, along=2)
copies = na.omit(copies[!is.na(rownames(copies)),])
cancer_copies = t(t(copies) / purity$estimate)
cancer_copies = cancer_copies[rownames(cancer_copies) %in% splice$ensembl_gene_id,]
clust_x = tibble(Sample = colnames(cancer_copies),
                 umap_x = rank(uwot::umap(t(cancer_copies), n_components=1)[,1]))
clust_y = tibble(ensembl_gene_id = rownames(cancer_copies),
                 umap_y = rank(uwot::umap(cancer_copies, n_components=1)[,1]))

# load compensation scores
comp = readRDS()

dset = reshape2::melt(cancer_copies, value.name="cancer_copies") %>%
    inner_join(clust_x) %>%
    inner_join(clust_y) %>%
    inner_join(splice) %>%
#    plt$cluster(cancer_copies ~ sid + Sample) %>%
    mutate(cancer_copies = pmin(cancer_copies, 4),
           hgnc_chr = sprintf("%s -- %s", hgnc_symbol, chr))

# umap plot for different samples & genes
pdf("do.pdf", 20, 12)
ggplot(dset, aes(x=forcats::fct_reorder(Sample, umap_x),
                 y=forcats::fct_reorder(hgnc_chr, umap_y), fill=cancer_copies)) +
    geom_tile() +
    scale_fill_distiller(palette="Spectral") +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size=5),
          strip.text.y.left = element_text(angle = 0)) +
    ggtitle("copy numbers")

ggplot(dset, aes(x=forcats::fct_reorder(Sample, umap_x),
                 y=forcats::fct_reorder(ensembl_gene_id, umap_y), fill=cancer_copies)) +
    geom_tile() +
    scale_fill_distiller(palette="Spectral") +
    facet_grid(subunit ~ ., space="free", scales="free", switch="y") +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          strip.text.y.left = element_text(angle = 0)) +
    ggtitle("copy numbers")
dev.off()
