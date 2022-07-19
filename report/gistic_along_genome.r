library(dplyr)
library(ggplot2)
seq = import('seq')
tcga = import('data/tcga')

hlg = c("MYC", "EGFR", "RBM14", "CDKN1A", "TP53")
gt = seq$gene_table() %>%
    transmute(chr = factor(chromosome_name, levels=c(1:22,'X')),
              gene_name = external_gene_name,
              tss = transcription_start_site) %>%
    filter(!is.na(chr)) %>%
    group_by(chr, gene_name) %>%
        summarize(tss = mean(tss)) %>%
    ungroup()
gt2 = gt %>% filter(gene_name %in% hlg)

gg = tcga$cna_gistic(thresh=TRUE)
cg = tibble(gene_name = rownames(gg),
            f_amp = rowSums(gg > 0) / nrow(gg),
            f_del = -rowSums(gg < 0) / nrow(gg)) %>%
    inner_join(gt) %>%
    filter(chr %in% c(1:22,'X')) %>%
    mutate(chr = factor(chr, levels=c(1:22,'X')))

cg2 = tidyr::gather(cg, "type", "frac", -chr, -gene_name, -tss)

ggplot(cg2, aes(x=tss)) +
    geom_hline(yintercept=0, color="black") +
    geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
#    geom_line(data=cg2 %>% filter(abs(frac)>0.1), aes(y=frac, group=type, color=type), size=2) +
    scale_fill_manual(values=c(f_amp="firebrick", f_del="navy")) +
    geom_vline(data=gt2, aes(xintercept=tss), linetype="dashed", color="grey") +
    geom_point(data=gt2, y=0) +
    geom_text(data=gt2, aes(label=paste(" ", gene_name)), y=0, hjust=0, vjust=-0.5, angle=45, size=3) +
    facet_grid(. ~ chr, scales="free", space="free") +
    theme_minimal()
#todo: high/low amp
