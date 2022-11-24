library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga')
cm = import('./common')

schema = function() {
    img = grid::rasterGrob(magick::image_read("external/schema.svg", density=300))
    ggplot() + annotation_custom(img) + theme(panel.background=element_blank())
}

overlap = function() {
    img = grid::rasterGrob(magick::image_read("external/overlap2.svg"))
    ggplot() + annotation_custom(img) + theme(panel.background=element_blank())
}

cna_along_genome = function(gistic) {
    labs = gistic$genes %>% filter(gene_name %in% cm$hlg) %>%
        group_by(gene_name) %>%
            slice_max(abs(frac)) %>%
        ungroup() %>%
        inner_join(gistic$smooth %>% select(type, chr, gam)) %>%
        rowwise() %>%
        mutate(frac = mgcv::predict.gam(gam, newdata=data.frame(tss=tss)))

    smooth = gistic$smooth %>% select(-gam) %>% tidyr::unnest(steps) %>%
        mutate(frac_amp = ifelse(frac > 0.15, frac, NA),
               type = stringr::str_to_title(type))
    ggplot(smooth, aes(x=tss)) +
        geom_hline(yintercept=0, color="black") +
        geom_hline(yintercept=0.15, color="black", linetype="dashed") +
        geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
        scale_fill_manual(values=cm$cols[c("Amplification", "Deletion")], name="CNA") +
        geom_line(aes(y=frac_amp, group=type, color="Frequently\namplified"),
                  lineend="round", size=1) +
        scale_color_manual(values=c("Frequently\namplified"="#960019"), name="") +
        geom_point(data=labs, aes(y=frac), color="black", fill="white", shape=21) +
        ggrepel::geom_text_repel(data=labs, aes(y=frac, label=gene_name), size=3,
                                 point.size=10, max.iter=1e5, max.time=10) +
        facet_grid(. ~ chr, scales="free", space="free") +
        labs(y = "Alteration frequency TCGA") +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.spacing.x = unit(1, "mm")) +
        coord_cartesian(clip="off", expand=FALSE)
}

cna_expr_scales = function() {
    dset = readRDS("../data/df_ccle.rds") %>%
        group_by(gene, covar) %>%
            filter(n() > 30, sum(expr == 0) < 0.1, mean(expr) > 10, mean(abs(copies-2) > 1) > 0.15) %>%
            mutate(expr = expr/sf, expr = expr/mean(expr[abs(2-copies)<0.15])) %>%
        ungroup()

    ggplot(dset, aes(x=copies, y=expr)) +
        geom_hex(aes(color=..count..), bins=50) +
        geom_hline(yintercept=0, color="white", linetype="dashed") +
        geom_vline(xintercept=0, color="white", linetype="dashed") +
        geom_density2d(color="black", bins=8, size=0.3) +
        scale_color_distiller(palette="Greys", direction=1, guide="none") +
        scale_fill_distiller(palette="Greys", direction=1, breaks=c(100,1000)) +
        xlim(c(0, 4)) + ylim(c(0,2)) +
        geom_abline(slope=0.5, intercept=0, color="blue", linetype="dashed") +
        labs(x="DNA copies", y="Normalized\nexpression CCLE") +
        theme_classic() +
        coord_fixed(ratio=2, expand=FALSE)
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")
    cosmic = cm$get_cosmic_annot()

    top = cna_along_genome(gistic)
    left = (wrap_elements(cna_expr_scales()) /
            wrap_elements(overlap() + theme(plot.margin=margin(0,0,0,-15,"mm")))) +
        plot_layout(heights=c(2.5,3))
    right = schema()

    asm = (top / ((left | right) + plot_layout(widths=c(1,2)))) +
        plot_layout(heights=c(1,3)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig1-Motivation.pdf", 14, 10)
    print(asm)
    dev.off()
})
