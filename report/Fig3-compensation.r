library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
fig1 = import('./Fig1-Motivation')

schema = function() {
    schema = grid::rasterGrob(magick::image_read("external/compensation.svg"))
    ggplot() + annotation_custom(schema) + theme(panel.background=element_blank())
}

tcga_ccle_cor = function(gistic_amp, cosmic) {
    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    both = inner_join(ccle, tcga3, by="gene") %>%
        dplyr::rename(gene_name = gene) %>%
        left_join(cosmic) %>%
        inner_join(gistic_amp)

    dx = ggplot(both, aes(x=estimate.x)) +
        geom_density(fill="#dedede") +
        theme_void() +
        scale_y_continuous(expand=c(0,0))
    dy = ggplot(both, aes(x=estimate.y)) +
        geom_density(fill="#dedede") +
        theme_void() +
        scale_x_continuous(expand=c(0,0)) +
        coord_flip(expand=FALSE) +
        plot_layout(tag_level="new")

    m = lm(estimate.y ~ estimate.x, data=both) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.2f~\n~p~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    p = ggplot(both, aes(x=estimate.x, y=estimate.y)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_rect(xmin=-Inf, ymin=-Inf, xmax=-0.3, ymax=-0.3, color="orange",
                  linetype="dashed", aes(fill="Compensated")) +
        scale_fill_manual(values=c(Compensated="#fff9f5de"), name="Area") +
        geom_point(data=both %>% filter(is.na(type)),
                   aes(size=n_aneup.y), alpha=0.2) +
        geom_density2d(color="white", breaks=c(0.5, 0.99)) +
        geom_point(data=both %>% filter(!is.na(type)),
                   aes(size=n_aneup.y, color=type), alpha=0.7) +
        geom_smooth(method="lm", color="blue") +
        annotate("text", y=1.3, x=-0.9, hjust=0, label=lab, color="blue", parse=TRUE) +
        ggrepel::geom_label_repel(data=both %>% filter(!is.na(type)),
                   aes(label=gene_name, color=type), max.overlaps=11,
                   size=3, min.segment.length=0,
                   segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        theme_classic() +
        labs(x = "Expression over expected amplified genes CCLE",
             y = "Expression over expected amplified genes TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
}

og_tsg_comp = function(gistic_amp, cosmic) {
    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))

    dset = list(CCLE=ccle, TCGA=tcga3) %>% bind_rows(.id="dset") %>%
        dplyr::rename(gene_name = gene) %>%
        left_join(cosmic) %>%
        inner_join(gistic_amp) %>%
        mutate(type = ifelse(is.na(type), "Background", type))
    meds = dset %>% group_by(dset) %>%
        summarize(bg = median(estimate[type == "Background"], na.rm=TRUE))

    ggplot(dset, aes(x=type, y=estimate, color=type)) +
        geom_boxplot(outlier.shape=NA) +
        ggbeeswarm::geom_quasirandom(alpha=0.5) +
        ggsignif::geom_signif(y_position=c(1.6, 1.8), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))) +
        geom_hline(data=meds, aes(yintercept=bg), linetype="dashed", color="black") +
        facet_wrap(~ dset) +
        theme_classic() +
        ggtitle("Compensation of Oncogenes/TSGs")
}

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
    gistic_amp = fig1$get_gistic_scores() %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)
    cosmic = fig1$get_cosmic_annot()
    all = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")

    left = (schema() / tcga_ccle_cor(gistic_amp, cosmic)) +
        plot_layout(heights=c(1,4))

    right = (og_tsg_comp(gistic_amp, cosmic) / comp_orf(all, gistic_amp)) +
        plot_layout(heights=c(1,2))

    asm = (left | right) + plot_layout(widths=c(4,3)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig3-compensation.pdf", 15, 10)
    print(asm)
    dev.off()
})
