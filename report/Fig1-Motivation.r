library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga')

get_gistic_scores = function() {
    gt = seq$gene_table() %>%
        transmute(chr = factor(chromosome_name, levels=c(1:22,'X')),
                  gene_name = external_gene_name,
                  tss = transcription_start_site) %>%
        filter(!is.na(chr)) %>%
        group_by(chr, gene_name) %>%
            summarize(tss = mean(tss)) %>%
        ungroup()

    readRDS("../data/tcga_prod2-gistic.rds") %>%
        transmute(gene_name = SYMBOL,
                  type = ALT_TYPE,
                  frac = ifelse(type == "amplification", OVERALL_FREQ, -OVERALL_FREQ)) %>%
        inner_join(gt)
}

cna_along_genome = function(gistic_scores, hlg=c()) {
    labs = gistic_scores %>% filter(gene_name %in% hlg) %>%
        group_by(gene_name) %>%
            slice_max(abs(frac)) %>%
        ungroup()

    ggplot(cg2, aes(x=tss)) +
        geom_hline(yintercept=0, color="black") +
        geom_rect(xmin=-Inf, xmax=Inf, ymin=-0.15, ymax=0.15, fill="#efefef", color=NA) +
        geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
        scale_fill_manual(values=c(amplification="firebrick", deletion="navy"), name="CNA") +
        geom_vline(data=labs, aes(xintercept=tss), linetype="dashed", color="grey") +
        ggrepel::geom_text_repel(data=labs, aes(y=frac, label=gene_name),
                                 min.segment.length=0) +
        facet_grid(. ~ chr, scales="free", space="free") +
        labs(y = "Alteration frequency") +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.spacing.x = unit(1, "mm")) +
        coord_cartesian(clip="off")
}

get_cosmic_annot = function() {
    manual = readRDS("../data/genesets/manual.rds")
    manual[grepl("Cosmic", names(manual))] %>%
        stack() %>% as_tibble() %>%
        mutate(ind = ifelse(grepl("OG_Hallmark", ind), paste(ind, "Tier1"), as.character(ind))) %>%
        transmute(gene_name = values,
                  type = case_when(
                      grepl("OG", ind) ~ "Oncogene",
                      grepl("TSG", ind) ~ "TSG"
                  ),
                  tier = sub(".*Tier([12])$", "\\1", ind)) %>%
        distinct()
}

og_vs_tsg = function(gistic, cosmic, hlg=c()) {
    cosmic2 = cosmic %>% group_by(gene_name, tier) %>%
        summarize(type = ifelse(length(type) == 1, type, "both"))
    gwide = tidyr::pivot_wider(gistic, names_from="type", values_from="frac") %>%
        left_join(cosmic2) %>%
        mutate(label = ifelse(gene_name %in% hlg, gene_name, NA))

    p = ggplot(gwide, aes(x=amplification, y=-deletion, color=type)) +
        geom_abline(intercept=0, slope=1, color="grey", linetype="dashed") +
        geom_point(aes(shape=tier), na.rm=TRUE) +
        ggrepel::geom_label_repel(aes(label=label), size=3, na.rm=TRUE) +
        scale_shape_manual(values=c("1"=19, "2"=1), name="Tier") +
        scale_color_manual(values=c(Oncogene="firebrick", TSG="navy", both="purple"), name="Type") +
        coord_fixed() +
        theme_classic() +
        labs(x = "Amplification frequency",
             y = "Deletion frequency")

    gwide$type[is.na(gwide$type)] = "background"
    boxbee = function(y1, y2) list(
        geom_boxplot(),
        ggbeeswarm::geom_quasirandom(alpha=0.5),
        ggsignif::geom_signif(comparisons=list(c("background", "Oncogene"), c("background", "TSG")),
                              y_position=c(y1, y2), color="black", test=wilcox.test),
        theme_classic(),
        theme(legend.position="none")
    )
    pa = ggplot(gwide, aes(x=type, y=amplification, color=type)) +
        boxbee(0.48, 0.52) + ylab("Amplification frequency")
    pd = ggplot(gwide, aes(x=type, y=-deletion, color=type)) +
        boxbee(0.4, 0.44) + ylab("Deletion frequency")

    (p | pa | pd) + plot_layout(widths=c(5,2,2))
}

sys$run({
    hlg = c("MYC", "EGFR", "CCND1", "CDKN1A", "TP53", "BAP1", "CDKN1A", "IL7R", "CKS1B",
            "APC", "CDKN2A", "KRAS", "NRAS", "RB1", "SMAD4", "CCNE1", "PIK3CA", "AURKA")

    gistic = get_gistic_scores()
    cosmic = get_cosmic_annot()

    top = cna_along_genome(gistic, hlg)
    btm = og_vs_tsg(gistic, cosmic, hlg)

    asm = (top / btm) + plot_layout(heights=c(3,4)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig1-Motivation.pdf", 14, 8)
    print(asm)
    dev.off()
})
