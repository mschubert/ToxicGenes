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

tcga_ccle_cor = function() {
    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    both = inner_join(ccle, tcga3, by="gene") %>%
        left_join(fig1$get_cosmic_annot() %>% dplyr::rename(gene=gene_name))

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

    p = ggplot(both, aes(x=estimate.x, y=estimate.y)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_rect(xmin=-Inf, ymin=-Inf, xmax=-0.3, ymax=-0.3, color="orange",
                  linetype="dashed", color="orange", aes(fill="Compensated")) +
        scale_fill_manual(values=c(Compensated="#ffdead35"), name="Area") +
        geom_point(data=both %>% filter(is.na(type)),
                   aes(size=n_aneup.y), alpha=0.2) +
        geom_density2d(color="white", breaks=c(0.5, 0.99)) +
        geom_point(data=both %>% filter(!is.na(type)),
                   aes(size=n_aneup.y, color=type), alpha=0.7) +
        ggrepel::geom_label_repel(data=both %>% filter(!is.na(type)),
                   aes(label=gene, color=type), max.overlaps=11,
                   size=3, min.segment.length=0,
                   segment.alpha=0.3, fill="#ffffffc0", label.size=NA) +
        theme_classic() +
        labs(x = "CCLE", y = "TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
}

sys$run({
    (schema() / tcga_ccle_cor()) +
        plot_layout(heights=c(1,4)) +
        plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=18, face="bold"))
})
