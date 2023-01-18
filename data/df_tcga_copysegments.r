library(dplyr)
library(plyranges)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')
seq = import('seq')

# this does not really work because small changes in segment identify new
# segment, e.g. 3.01 to 3.05 copies

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'df_tcga_copysegments.rds')
)

cohorts = tcga$cohorts()
# excl x,y chroms

genes = seq$gene_table() %>%
    filter(gene_biotype == "protein_coding") %>%
    group_by(ensembl_gene_id) %>%
    summarize(seqnames=unique(chromosome_name), start=min(start_position),
        stop=max(end_position), strand=c("-","+")[unique(strand/2+1.5)]) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)

segs = readr::read_tsv("./gistic/TCGA.all_cancers.150601.zigg_events.160923.txt") %>%
    transmute(Sample=sample, seqnames=chr, start=base_start, end=base_end,
              seg_id = seq_len(nrow(.))) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)

res = join_overlap_intersect(segs, genes) %>%
    group_by(seg_id) %>%
    summarize(n_genes = n_distinct(ensembl_gene_id)) %>%
    as.data.frame() %>% as_tibble()

saveRDS(res, file=args$outfile)
