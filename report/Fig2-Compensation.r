library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

schema = function() {
    img = grid::rasterGrob(magick::image_read("external/comp.svg"))
    ggplot() + annotation_custom(img) + theme(panel.background=element_blank())
}

tcga_ccle_cor = function(both, gistic_amp, cosmic) {
    dx = ggplot(both, aes(x=estimate.x)) +
        geom_density(fill=cm$cols["CCLE"], alpha=0.5) +
        theme_void() +
        scale_y_continuous(expand=c(0,0))
    dy = ggplot(both, aes(x=estimate.y)) +
        geom_density(fill=cm$cols["TCGA"], alpha=0.5) +
        theme_void() +
        scale_x_continuous(expand=c(0,0)) +
        coord_flip(expand=FALSE) +
        plot_layout(tag_level="new")

    both$group = with(both, case_when(
        estimate.x < -0.3 & estimate.y < -0.3 ~ "Compensated",
        estimate.x > 0.3 & estimate.y > 0.3 ~ "Hyperactivated",
        TRUE ~ "Background"
    ))
    both$pdist = with(both, sqrt(estimate.x^2 + estimate.y^2))

    hl = c("RBM14", "CDKN1A", "SNRPA", "ZBTB14", "POU2F1", "STAT3", "BUB1",
           "BRCA1", "RBM33", "DNMT3A", "CCNE1", "MDM2", "ERBB2", "CDC73",
           "SRSF3", "AKT3", "BAX", "POLQ", "MPP2", "CYP1B1-AS1", "ZNF879",
           "ID2", "ZNF418")
    cols = cm$cols[c("Compensated", "Hyperactivated", "Background")]
    m = lm(estimate.y ~ estimate.x, data=both) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.2f~\n~p~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    p = ggplot(both, aes(x=estimate.x, y=estimate.y)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_point(data=both %>% filter(group == "Background"),
                   aes(size=n_aneup.y, color=group), alpha=0.2) +
        geom_point(data=both %>% filter(group != "Background"),
                   aes(size=n_aneup.y, color=group), alpha=0.7) +
        scale_color_manual(values=cols, name="Compensation\nstatus") +
        scale_size_area(max_size=6, breaks=c(100, 500, 1000)) +
        annotate("text", x=0.05, y=-0.8, label="Compensated", color=cm$cols[["Compensated"]],
                 size=6, fontface="bold", alpha=0.7, hjust=0) +
        annotate("text", x=0.4, y=1.55, label="Hyperactivated", color=cm$cols[["Hyperactivated"]],
                 size=4, fontface="bold", alpha=0.7, hjust=0) +
        ggrepel::geom_label_repel(data=both %>% filter(pdist > 1 | gene_name %in% hl),
                   aes(label=gene_name, color=group), max.overlaps=20,
                   box.padding=unit(0.1, "lines"),
                   size=3, min.segment.length=0, max.iter=1e6, max.time=30,
                   segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        theme_classic() +
        labs(size = "TCGA\nAmplifications",
             x = "Δ Expression over expected CCLE",
             y = "Δ Expression over expected TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
}

go_tcga_ccle = function(n=5) {
    ccle = readRDS("../ccle/pan/stan-nb/all/GO_Biological_Process_2021.rds")$amp %>%
        transmute(label=label, stat_CCLE=statistic, fdr_CCLE=adj.p)
    tcga = readRDS("../tcga/pan/stan-nb_puradj/all/GO_Biological_Process_2021.rds")[[1]] %>%
        transmute(label=label, stat_TCGA=statistic, fdr_TCGA=adj.p)
    dset = inner_join(ccle, tcga) %>%
        mutate(mean_stat = rowMeans(cbind(stat_CCLE, stat_TCGA))) %>%
        arrange(mean_stat)
    dset = bind_rows(Compensated=head(dset, n), Hyperactivated=tail(dset,n), .id="type") %>%
        tidyr::pivot_longer(c(stat_CCLE, stat_TCGA),
                            names_pattern="_(TCGA|CCLE)",
                            names_to="dset")

    plot_one = function(dset, .desc) {
        dset$label = forcats::fct_reorder(dset$label, dset$mean_stat, .desc=.desc)
        txt = dset %>% select(label) %>% distinct()
        ggplot(dset, aes(x=label, y=value)) +
            geom_col(aes(fill=dset), width=1.2, position=position_dodge(width=0.6), alpha=0.5) +
            geom_text(data=txt, aes(label=paste(" ", label)), y=0, hjust=0) +
            coord_flip(expand=FALSE, clip="off") +
            scale_fill_manual(values=cm$cols[c("TCGA", "CCLE")], name="Dataset") +
            labs(x = "Gene Ontology", y = "Δ Expression over expected in CCLE/TCGA (Wald statistic)") +
            theme_classic() +
            theme(axis.text.y = element_blank())
    }
    p1 = plot_one(dset %>% filter(type == "Compensated"), TRUE) +
        ggtitle("Compensated") + scale_y_reverse() + theme(axis.title.x = element_blank())
    p2 = plot_one(dset %>% filter(type == "Hyperactivated"), FALSE) +
        ggtitle("Hyperactivated") + plot_layout(tag_level="new")
    (p1 / p2) + plot_layout(guides="collect")
}

cna_comp = function(gistic, comp_all) {
    gwide = gistic %>%
        tidyr::pivot_wider(names_from="type", values_from="frac")
    gwide = gwide %>%
        mutate(type = case_when(
            amplification > 0.15 & deletion < -0.15 ~ "Amp+Del",
            amplification > 0.15 ~ "Amplified",
            deletion < -0.15 ~ "Deleted",
            TRUE ~ NA_character_
        )) %>% filter(!is.na(type)) %>% bind_rows(gwide %>% mutate(type="Background"))

    both = inner_join(comp_all %>% select(-type), gwide) %>%
        mutate(type = factor(type, levels=c("Background", "Amplified", "Deleted"))) %>%
        filter(!is.na(type))

    common = function(y, coordy, sigy) list(
        geom_boxplot(outlier.shape=NA, alpha=0.7),
        ggsignif::geom_signif(y_position=sigy, color="black", test=t.test,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))),
        scale_fill_manual(values=cm$cols[c("Background", "Amplified", "Deleted")]),
        labs(fill = "Frequent CNA", x = "Copy number subset", y = "Δ ORF (Wald statistic)"),
        theme_classic(),
        coord_cartesian(ylim=coordy),
        theme(axis.text.x = element_blank()),
        geom_hline(yintercept=median(y[both$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
    )

    p1 = ggplot(both, aes(x=type, y=estimate.x, fill=type)) +
        common(both$estimate.x, c(-0.3, 0.7), c(0.4, 0.55)) +
        labs(title = "CCLE", y="Δ Expression / expected")
    p2 = ggplot(both, aes(x=type, y=estimate.y, fill=type)) +
        common(both$estimate.y, c(-0.5, 1.4), c(0.95, 1.2)) +
        labs(title = "TCGA", y="")

    (p1 | (p2 + plot_layout(tag_level="new"))) + plot_layout(guides="collect")
}

og_comp = function(comp) {
    both = comp %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG"))) %>%
        filter(!is.na(type))

    common = function(y, coordy, sigy) list(
        geom_boxplot(outlier.shape=NA, alpha=0.7),
        ggsignif::geom_signif(y_position=sigy, color="black", test=t.test,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))),
        scale_fill_manual(values=cm$cols[c("Background", "Oncogene", "TSG")]),
        labs(fill = "Driver status\n(freq. amplified)", x = "Gene type subset"),
        theme_classic(),
        coord_cartesian(ylim=coordy),
        theme(axis.text.x = element_blank()),
        geom_hline(yintercept=median(y[both$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
    )

    p1 = ggplot(both, aes(x=type, y=estimate.x, fill=type)) +
        common(both$estimate.x, c(-0.3, 0.7), c(0.4, 0.55)) +
        labs(title = "CCLE", y="Δ Expression / expected")
    p2 = ggplot(both, aes(x=type, y=estimate.y, fill=type)) +
        common(both$estimate.y, c(-0.5, 1.4), c(1.0, 1.2)) +
        labs(title = "TCGA", y="")

    (p1 | (p2 + plot_layout(tag_level="new"))) + plot_layout(guides="collect")
}

sys$run({
    cosmic = cm$get_cosmic_annot()
    gistic = readRDS("../data/gistic_smooth.rds")$genes
    gistic_amp = gistic %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)

    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    comp_all = inner_join(ccle, tcga3, by="gene") %>%
        dplyr::rename(gene_name = gene) %>%
        left_join(cosmic)
    comp = comp_all %>% inner_join(gistic_amp)

    left = (wrap_elements(schema() + theme(plot.margin=margin(0,5,0,-10,"mm")))) /
        tcga_ccle_cor(comp, gistic_amp, cosmic)
    right = cna_comp(gistic, comp_all) / og_comp(comp) / go_tcga_ccle()

    asm = ((left + plot_layout(heights=c(1,3))) |
        (right + plot_layout(heights=c(1,1,2)))) + plot_layout(widths=c(3.5,2)) +
        plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("Fig2-Compensation.pdf", 14, 10)
    print(asm)
    dev.off()
})
