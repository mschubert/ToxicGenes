library(dplyr)
library(ggplot2)
seq = import('seq')
tcga = import('data/tcga')

#hlg = c("MYC", "EGFR", "RBM14", "CDKN1A", "TP53")
hlg = readRDS("../data/genesets/manual.rds")[c("Davoli_oncogenes", "Davoli_TSGs")] %>%
    unlist(use.names=FALSE) %>% na.omit() %>% c()

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
    geom_point(data=gt2, aes(color=group), y=0) +
    ggrepel::geom_label_repel(data=gt2, aes(color=group,label=gene_name), y=0, hjust=0, vjust=-0.5, angle=0, size=5, max.time=10, max.iter=1e5) +
    facet_grid(. ~ chr, scales="free", space="free") +
    theme_minimal() +
    coord_cartesian(clip="off")
#todo: high/low amp


pos_set_comp = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv") %>%
    dplyr::rename(gene_name=gene)
hit = pos_set_comp %>% filter(hit)
# comp selection
gt2 = inner_join(hit, cg)

# davoli
dav = readRDS("../data/genesets/manual.rds")[c("Davoli_oncogenes", "Davoli_TSGs")] %>%
    lapply(function(x) x[!is.na(x)]) %>% stack() %>% dplyr::rename(gene_name=values, group=ind)
gt2 = inner_join(dav, cg) %>% as_tibble()
ggplot(gt2, aes(x=f_amp, y=f_del)) +
    geom_point(aes(color=group)) +
    ggrepel::geom_text_repel(aes(color=group,label=gene_name), size=3, max.time=10, max.iter=1e5)

xx = right_join(dav, pos_set_comp) %>% inner_join(cg) %>% as_tibble() %>%
    select(gene_name, group, est_ccle, est_tcga, f_amp, f_del)
xxl = xx %>%
    tidyr::gather("dset", "estimate", -gene_name, -group, -f_amp, -f_del)
p1 = ggplot(xxl, aes(x=dset, fill=group, y=estimate)) +
    geom_boxplot() + ggtitle("all")
p2 = ggplot(xxl %>% filter(f_amp > 0.05), aes(x=dset, fill=group, y=estimate)) +
    geom_boxplot() + ggtitle("min 5% amp")
p3 = ggplot(xxl %>% filter(f_del < -0.05), aes(x=dset, fill=group, y=estimate)) +
    geom_boxplot() + ggtitle("min 5% del")
library(patchwork)
(p1 | p2 | p3) + plot_layout(guides="collect")


yy = xx %>% filter(f_amp > 0.05) %>% mutate(label=ifelse(gene_name %in% hit$gene_name, gene_name, NA))
ggplot(yy, aes(x=est_ccle, y=est_tcga)) +
    geom_point() +
    geom_smooth(method="lm") +
    geom_density2d(color="green", bins=20)
#    ggrepel::geom_text_repel(aes(label=label))
