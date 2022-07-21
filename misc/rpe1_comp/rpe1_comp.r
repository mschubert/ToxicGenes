library(dplyr)
library(DESeq2)
library(ggplot2)
seq = import('seq')

gt = seq$gene_table() %>%
    group_by(chr=chromosome_name, ensembl_gene_id, gene=external_gene_name) %>%
        summarize(start = mean(start_position)) %>%
    ungroup() %>%
    filter(chr %in% c(1:22,'X'))

comp = readr::read_tsv("../../cor_tcga_ccle/positive_comp_set.tsv") %>%
    filter(hit) %>% select(-hit) %>%
    left_join(gt) %>%
    arrange(gene) %>%
    select(ensembl_gene_id, chr, start, gene, everything())
write.table(comp, sep="\t", row.names=FALSE, quote=FALSE, file="comp_genes.tsv")

rpe = readRDS("../../../dorine/data/rnaseq.rds")
nm = rpe$normMatrix

eset = DESeq2::estimateSizeFactors(rpe$eset, normMatrix=nm) %>%
    counts(normalized=TRUE)
narray::intersect(eset, comp$ensembl_gene_id, nm, along=1)
names(dimnames(eset)) = c("ensembl_gene_id", "short")

dset = inner_join(reshape2::melt(eset, value.name="reads"),
                  reshape2::melt(nm, value.name="copies")) %>%
    inner_join(comp) %>% as_tibble() %>%
    filter(copies %in% 1:3) %>%
    group_by(gene) %>%
        filter(sum(copies %in% c(1,3)) > 0) %>%
    ungroup()

diffs = dset %>% group_by(gene) %>%

ggplot(dset, aes(x=gene, y=reads, color=factor(copies))) +
    ggbeeswarm::geom_quasirandom(data=dset%>% filter(copies==2), alpha=0.7) +
    ggbeeswarm::geom_quasirandom(data=dset%>% filter(copies!=2), alpha=0.7) +
    theme(axis.text.x = element_text(angle=60, hjust=1)) +
    scale_y_continuous(trans="log1p", breaks=c(0,1,2,5,10,50,200,1000,5000,20000)) +
    scale_color_manual(values=c("1"="navy", "2"="grey", "3"="firebrick"))
