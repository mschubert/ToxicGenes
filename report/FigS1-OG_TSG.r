library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
cm = import('./common')

og_vs_tsg = function(gistic, cosmic) {
    gwide = tidyr::pivot_wider(gistic, names_from="type", values_from="frac") %>%
        left_join(cosmic) %>%
        mutate(label = ifelse(!is.na(type) & gene_name %in% cm$hlg, gene_name, NA))

    areas = data.frame(xmin=c(0.15,-Inf), xmax=c(Inf,Inf), ymin=c(-Inf,0.15), ymax=c(Inf,Inf),
        status=c("Deleted", "Amplified"))

    ggplot(gwide, aes(x=-deletion, y=amplification, color=type)) +
        geom_rect(data=areas, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=status),
            linetype="dashed", alpha=0.1, inherit.aes=FALSE) +
        scale_fill_manual(values=cm$cols[c("Amplified", "Deleted")], name="Frequently") +
        geom_point(aes(shape=tier), na.rm=TRUE) +
        ggrepel::geom_label_repel(aes(label=label), size=3, min.segment.length=0,
            segment.alpha=0.3, fill="#ffffffc0", label.size=NA, na.rm=TRUE) +
        scale_shape_manual(values=c("1"=19, "2"=1), name="Tier") +
        scale_color_manual(values=cm$cols[c("Oncogene", "TSG", "OG+TSG")], name="Driver status") +
        coord_fixed() +
        theme_classic() +
        labs(x = "Deletion frequency",
             y = "Amplification frequency")
}

venn_amp_del = function(gistic, cosmic) {
    amps = gistic %>% filter(type == "amplification")
    dels =  gistic %>% filter(type == "deletion")
    sets = list(
        `All genes` = unique(gistic$gene_name),
        `Frequently\namplified` = unique(amps$gene_name[amps$frac > 0.15]),
        `Frequently\ndeleted` = unique(dels$gene_name[dels$frac < -0.15])
    )
    cols = cm$cols[c("Background", "Amplified", "Deleted")]
    names(cols) = names(sets)
    plt$venn(sets, alpha=0.3) +
        scale_fill_manual(values=cols)
}

venn_og_tsg = function(gistic, cosmic) {
    amps = gistic %>% filter(type == "amplification")
    dels =  gistic %>% filter(type == "deletion")
    sets = list(
        `Frequently\namplified` = unique(amps$gene_name[amps$frac > 0.15]),
        Oncogene = cosmic$gene_name[grepl("^O", cosmic$type)],
        TSG = cosmic$gene_name[grepl("TSG", cosmic$type)]
    )
    cols = cm$cols[c("Oncogene", "TSG", "Amplified")]
    names(cols)[3] = names(sets)[1]
    plt$venn(sets, alpha=0.3) +
        scale_fill_manual(values=cols)
}

# DE normal-cancer for OG/TSG

# RPE-1 common parental DE with amps?

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")$genes
    cosmic = cm$get_cosmic_annot()

    btm = (og_vs_tsg(gistic, cosmic) | venn_amp_del(gistic, cosmic) | venn_og_tsg(gistic, cosmic)) +
        plot_layout(widths=c(3,2,2))

    asm = btm + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS1-OG_TSG.pdf", 14, 5.5)
    print(asm)
    dev.off()
})
