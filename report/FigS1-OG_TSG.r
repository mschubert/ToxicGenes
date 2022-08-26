library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga')
cm = import('./common')

og_vs_tsg = function(gistic, cosmic, hlg=c()) {
    gwide = tidyr::pivot_wider(gistic, names_from="type", values_from="frac") %>%
        left_join(cosmic) %>%
        mutate(label = ifelse(!is.na(type) & gene_name %in% hlg, gene_name, NA))

    p = ggplot(gwide, aes(x=amplification, y=-deletion, color=type)) +
        geom_abline(intercept=0, slope=1, color="grey", linetype="dashed") +
        geom_point(aes(shape=tier), na.rm=TRUE) +
        ggrepel::geom_label_repel(aes(label=label), size=3, min.segment.length=0,
            segment.alpha=0.3, fill="#ffffffc0", label.size=NA, na.rm=TRUE) +
        scale_shape_manual(values=c("1"=19, "2"=1), name="Tier") +
        scale_color_manual(values=cm$cols[c("Oncogene", "TSG", "OG+TSG")], name="Driver status") +
        coord_fixed() +
        theme_classic() +
        labs(x = "Amplification frequency",
             y = "Deletion frequency")

    gwide$type[is.na(gwide$type)] = "Background"
    gwide$type = factor(gwide$type, levels=c("Background", "Oncogene", "TSG", "OG+TSG"))
    boxbee = function(y1, y2) list(
        geom_boxplot(outlier.shape=NA),
        scale_fill_manual(values=cm$cols[levels(gwide$type)]),
        ggsignif::geom_signif(comparisons=list(c("Background", "Oncogene"), c("Background", "TSG")),
                              y_position=c(y1, y2), color="black", test=wilcox.test),
        coord_cartesian(ylim=c(NA, max(c(y1, y2))+0.05)),
        theme_classic(),
        labs(x = "Gene type"),
        theme(axis.text.x = element_blank())
    )
    pa = ggplot(gwide, aes(x=type, y=amplification, fill=type)) +
        boxbee(0.45, 0.52) + ylab("Frequency") + ggtitle("Amplifications") +
        geom_hline(yintercept=median(gwide$amplification[gwide$type=="Background"]),
                   linetype="dashed", color="black")
    pd = ggplot(gwide, aes(x=type, y=-deletion, fill=type)) +
        boxbee(0.4, 0.47) + ylab("Frequency") + ggtitle("Deletions") +
        geom_hline(yintercept=-median(gwide$deletion[gwide$type=="Background"]),
                   linetype="dashed", color="black")

    (p | (pa / pd) + plot_layout(guides="collect")) + plot_layout(widths=c(2,1))
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")$genes
    cosmic = cm$get_cosmic_annot()

    btm = og_vs_tsg(gistic, cosmic, cm$hlg)

    asm = btm + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS1-OG_TSG.pdf", 10, 5.5)
    print(asm)
    dev.off()
})
