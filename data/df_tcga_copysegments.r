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

purity = tcga$purity() %>%
    filter(!is.na(estimate)) %>%
    transmute(Sample = Sample,
              purity = estimate)

genes = seq$gene_table() %>%
    filter(gene_biotype == "protein_coding") %>%
    group_by(ensembl_gene_id) %>%
    summarize(seqnames=unique(chromosome_name), start=min(start_position),
        stop=max(end_position), strand=c("-","+")[unique(strand/2+1.5)]) %>%
    makeGRangesFromDataFrame(keep.extra.columns=TRUE)

segs = do.call(c, lapply(cohorts, tcga$cna_segments, granges=TRUE)) %>%
    mutate(seg_id = seq_len(.))
segs$purity = purity$purity[match(segs$Sample, purity$Sample)]
segs = segs[with(segs, !is.na(purity) & abs((ploidy/2-1)/purity) > 0.75)]
segs$seg_id = seq_along(segs)

res = join_overlap_intersect_within(segs, genes) %>%
    group_by(seg_id) %>%
    summarize(n_genes = n_distinct(ensembl_gene_id))

saveRDS(res, file=args$outfile)
