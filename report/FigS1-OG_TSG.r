library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
cm = import('./common')

og_tsg_cna = function(gistic, cosmic) {
    dset = gistic %>%
        mutate(frac = abs(frac),
               cna = stringr::str_to_title(type)) %>% select(-type) %>%
        left_join(cosmic) %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG"))) %>%
        filter(!is.na(type))

    bg_line = dset %>% group_by(cna) %>%
        summarize(frac = median(frac[type == "Background"], na.rm=TRUE))

    ggplot(dset, aes(x=type, y=frac, fill=type)) +
        geom_boxplot(outlier.shape=NA) +
        scale_fill_manual(values=cm$cols[levels(dset$type)]) +
        labs(y="Frequency TCGA", x ="Gene type subset",
             fill="Driver status\n(whole genome)") +
        geom_hline(data=bg_line, aes(yintercept=frac), linetype="dashed", color="black") +
        ggsignif::geom_signif(comparisons=list(c("Background", "Oncogene"), c("Background", "TSG")),
            y_position=c(0.43,0.48,0.43,0.48), color="black", test=t.test, textsize=3) +
        facet_wrap(~ cna) +
        coord_cartesian(ylim=c(0.02, 0.52)) +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              strip.text.x = element_text(size=12),
              axis.text.x = element_blank())
}

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
        `Frequently amplified` = unique(amps$gene_name[amps$frac > 0.15]),
        `Frequently deleted` = unique(dels$gene_name[dels$frac < -0.15])
    )
    sets$Neither = setdiff(unique(gistic$gene_name), unlist(sets, use.names=FALSE))
    cols = cm$cols[c("Amplified", "Deleted", "Background")]
    names(cols) = names(sets)
    plt$venn(sets, alpha=0.3) +
        scale_fill_manual(values=cols) +
        theme(plot.margin = margin(0,5,0,5, unit="mm"))
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
        scale_fill_manual(values=cols) +
        theme(plot.margin = margin(0,5,0,5, unit="mm"))
}

de_og_tsg = function() {
# DE normal-cancer for OG/TSG
}

# RPE-1 common parental DE with amps?

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")$genes
    cosmic = cm$get_cosmic_annot()

    left = (og_vs_tsg(gistic, cosmic) / og_tsg_cna(gistic, cosmic)) +
        plot_layout(heights=c(2,1))

    right = (venn_amp_del(gistic, cosmic) / venn_og_tsg(gistic, cosmic))
        plot_layout(heights=c(3,4))

    asm = (left | right) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS1-OG_TSG.pdf", 11, 8)
    print(asm)
    dev.off()
})
