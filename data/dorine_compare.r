library(dplyr)
sys = import('sys')
seq = import('seq')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'dorine_compare.rds')
)

segs = readRDS("../../dorine/data/dna_vs_rna-compensation.rds") %>%
    filter(clone %in% c("14.10", "14.16", "14.21"), condition == "late") %>%
    tidyr::unnest(combined) %>%
    group_by(clone, condition, seqnames, seg_id) %>%
        summarize(lfc_DNA = unique(log2FoldChange.x),
                  lfc_RNA = mean(log2FoldChange.y)) %>%
    ungroup() %>%
    tidyr::pivot_longer(c(lfc_DNA, lfc_RNA), names_to="type",
                        names_prefix="lfc_", values_to="lfc")
gt = seq$gene_table() %>%
    group_by(ensembl_gene_id, seqnames=chromosome_name) %>%
        summarize(loc = mean(c(start_position, end_position))) %>%
    ungroup()
diff_expr = readRDS("../../dorine/basic_diff_expr/clones_vs_parental-copynaive.rds") %>%
    filter(clone %in% c("14.10", "14.16", "14.21"), condition == "late") %>%
    select(clone, genes) %>%
    tidyr::unnest(genes) %>%
    inner_join(gt) %>%
    filter(lfcSE < 4, seqnames != "10") %>%
    mutate(seqnames = droplevels(factor(seqnames, levels=c(1:22,'X'))))

saveRDS(list(segs=segs, diff_expr=diff_expr), file=args$outfile)
