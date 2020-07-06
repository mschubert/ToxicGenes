library(dplyr)
library(plyranges)
library(ggplot2)
io = import('io')
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('i', 'infile', 'file to read from', 'TableS2D.xlsx'),
    opt('o', 'outfile', 'file to save to', 'focal.RData'))

genes_hg19 = seq$coords$gene(idtype="hgnc_symbol", assembly="GRCh37", granges=TRUE)

segs = readxl::read_xlsx(args$infile, skip=19) %>%
    setNames(make.names(colnames(.))) %>%
    select(racs=Identifier, cohort=Cancer.Type, chr, start, stop,
           cna=Recurrent..Amplification...Deletion, Contained.genes)

pan = segs %>%
    filter(cohort == "PANCAN", cna == "Amplification") %>%
    select(-cohort, -Contained.genes) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)

racs_genes = join_overlap_intersect(genes_hg19, pan)

#tcga_cn = ...

tcga_comp = readxl::read_xlsx("../../tcga/pan/rlm3_puradj.xlsx", "amp") %>%
    filter(name %in% racs_genes$external_gene_name) %>%
    arrange(statistic)

orf_comp = readxl::read_xlsx("../../orf/fits_naive.xlsx", "pan") %>%
    filter(`GENE SYMBOL` %in% racs_genes$external_gene_name) %>%
    arrange(statistic)

both = inner_join(
        tcga_comp %>% dplyr::rename(gene=name, tcga_comp=statistic),
        orf_comp %>% transmute(gene=`GENE SYMBOL`, orf_stat=statistic)) %>%
    filter(tcga_comp < 0, orf_stat < 0, rsq > 0.02) %>%
    mutate(tcga_comp = pmax(tcga_comp, -50),
           orf_stat = pmax(orf_stat, -10),
           drank = rank(tcga_comp/5 + orf_stat),
           label = ifelse(drank < 30, gene, NA))

both %>%
    filter(!is.na(label), rsq>0.05, orf_stat< -2) %>%
    arrange(gene) %>%
    as.data.frame()

ggplot(both, aes(x=tcga_comp, y=orf_stat)) +
    geom_point(aes(alpha=rsq)) +
    ggrepel::geom_text_repel(aes(label=label))
dev.off()
