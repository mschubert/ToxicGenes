library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
orf = import('../orf/fits')
fig1 = import('./Fig1-Motivation')

schema = function() {
    schema = grid::rasterGrob(magick::image_read("external/ORFscreen.svg"))
    p = ggplot() + annotation_custom(schema) + theme(panel.background=element_blank())
    wrap_elements(p)
}

orf_volc = function(orfdata) {
    plt$volcano(orfdata, label_top=35, pos_label_bias=3, max.overlaps=100)
}

og_tsg_orf = function(orfdata) {
    cosmic = fig1$get_cosmic_annot()
    both = left_join(orfdata, cosmic) %>%
        mutate(type = ifelse(is.na(type), "Background", type))

    ggplot(both, aes(x=type, y=statistic, color=type)) +
        geom_boxplot(outlier.shape=NA) +
        ggbeeswarm::geom_quasirandom() +
        ggsignif::geom_signif(y_position=c(6.5, 9), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))) +
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
        ggbeeswarm::geom_quasirandom(alpha=0.5) +
        ggsignif::geom_signif(y_position=c(6.5, 9), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

sys$run({
    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
        filter(gene_name != "LOC254896") # not in tcga/ccle data

    asm = (schema() / (((orf_volc(orfdata)) |
        (og_tsg_orf(orfdata) / amp_del_orf(orfdata))) + plot_layout(widths=c(3,2)))) +
            plot_layout(heights=c(1,2)) + plot_annotation(tag_levels='a') &
            theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig2-ORFscreen.pdf", 12, 10)
    print(asm)
    dev.off()
})
