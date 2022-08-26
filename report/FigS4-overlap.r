library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

comp_orf = function(all, gistic_amp) {
    dset = inner_join(all %>% dplyr::rename(gene_name=gene), gistic_amp) %>%
        mutate(est_ccle_tcga = (est_ccle + est_tcga)/2,
               dropout = stat_orf < -5,
               label = ifelse(hit & stat_orf < -5 | stat_orf < -12, gene_name, NA))

    m = lm(stat_orf ~ est_ccle_tcga, data=dset) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.3f~\n~p~`=`~%.2g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    ggplot(dset, aes(x=(est_ccle+est_tcga)/2, y=stat_orf)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_point(aes(color=hit, alpha=dropout)) +
        ggrepel::geom_label_repel(aes(label=label, color=hit), size=3,
            min.segment.length=0, segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        scale_alpha_manual(values=c("TRUE"=0.9, "FALSE"=0.2), name="Dropout") +
        annotate("text", y=8, x=0.7, hjust=0, label=lab, color="blue", parse=TRUE) +
        geom_smooth(method="lm") +
        theme_classic()
}

sys$run({
    gistic_amp = readRDS("../data/gistic_smooth.rds")$genes %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)
    cosmic = cm$get_cosmic_annot()
    all = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")

    asm = comp_orf(all, gistic_amp)

#    left = (schema() / tcga_ccle_cor(gistic_amp, cosmic)) +
#        plot_layout(heights=c(1,4))
#
#    right = (og_tsg_comp(gistic_amp, cosmic) / comp_orf(all, gistic_amp)) +
#        plot_layout(heights=c(1,2))
#
#    asm = (left | right) + plot_layout(widths=c(4,3)) + plot_annotation(tag_levels='a') &
#        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS4-overlap.pdf", 8, 6)
    print(asm)
    dev.off()
})
