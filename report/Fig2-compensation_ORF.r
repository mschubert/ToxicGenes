library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
fig1 = import('./Fig1-Motivation')

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
        labs(x = "Expression over expected CCLE",
             y = "Expression over expected TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
}

orf_volc = function(orfdata) {
    plt$volcano(orfdata, label_top=35, pos_label_bias=3, max.overlaps=100) +
        labs(size = "Number of ORFs")
}

og_tsg_orf = function(orfdata) {
    cosmic = fig1$get_cosmic_annot()
    both = left_join(orfdata, cosmic) %>%
        mutate(type = ifelse(is.na(type), "Background", type))

    ggplot(both, aes(x=type, y=statistic, color=type)) +
        geom_boxplot(outlier.shape=NA) +
        ggsignif::geom_signif(y_position=c(6.5, 9), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))) +
        coord_cartesian(ylim=c(-8, 11)) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

amp_del_orf = function(orfdata) {
    gwide = fig1$get_gistic_scores() %>%
        tidyr::pivot_wider(names_from="type", values_from="frac") %>%
        mutate(type = case_when(
            amplification > 0.15 & deletion < -0.15 ~ "Both",
            amplification > 0.15 ~ "Amplified",
            deletion < -0.15 ~ "Deleted",
            TRUE ~ "Background"
        ))

    both = inner_join(orfdata, gwide) %>%
        mutate(type = factor(type, levels=c("Background", "Both", "Amplified", "Deleted")))

    ggplot(both, aes(x=type, y=statistic, color=type)) +
        geom_boxplot(outlier.shape=NA) +
        ggsignif::geom_signif(y_position=c(5, 6.5), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))) +
        coord_cartesian(ylim=c(-5, 8)) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

sys$run({
    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
        filter(gene_name != "LOC254896") # not in tcga/ccle data

    gistic_amp = fig1$get_gistic_scores() %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)
    cosmic = fig1$get_cosmic_annot()

    orf_cors = (og_tsg_orf(orfdata) / amp_del_orf(orfdata)) & theme_classic()

    top = tcga_ccle_cor(gistic_amp, cosmic)
    btm = orf_volc(orfdata) | orf_cors

    asm = (top / btm) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig2-compensation_ORF.pdf", 8, 14)
    print(asm)
    dev.off()
})
