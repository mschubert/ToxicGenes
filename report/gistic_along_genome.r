library(dplyr)
library(ggplot2)
seq = import('seq')
tcga = import('data/tcga')

hlg = c("MYC", "EGFR", "CCND1", "CDKN1A", "TP53", "BAP1", "CDKN1A", "IL7R", "CKS1B",
        "APC", "CDKN2A", "KRAS", "NRAS", "RB1", "SMAD4", "CCNE1", "PIK3CA", "AURKA")
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

cg2 = readRDS("../data/tcga_prod2-gistic.rds") %>%
    transmute(gene_name = SYMBOL,
              type = ALT_TYPE,
              frac = ifelse(type == "amplification", OVERALL_FREQ, -OVERALL_FREQ)) %>%
    inner_join(gt)

gt2 = cg2 %>% filter(gene_name %in% hlg) %>%
    group_by(gene_name) %>%
        slice_max(abs(frac)) %>%
    ungroup()

ggplot(cg2, aes(x=tss)) +
    geom_hline(yintercept=0, color="black") +
    geom_rect(xmin=-Inf, xmax=Inf, ymin=-0.15, ymax=0.15, fill="#efefef", color=NA) +
    geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
    scale_fill_manual(values=c(amplification="firebrick", deletion="navy")) +
    geom_vline(data=gt2, aes(xintercept=tss), linetype="dashed", color="grey") +
    ggrepel::geom_text_repel(data=gt2, aes(y=frac, label=gene_name), size=5, min.segment.length=0) +
    facet_grid(. ~ chr, scales="free", space="free") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.spacing.x = unit(1, "mm")) +
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
