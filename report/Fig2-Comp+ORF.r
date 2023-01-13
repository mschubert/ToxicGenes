library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

schema_comp = function() {
    img = grid::rasterGrob(magick::image_read("external/comp.svg"))
    ggplot() + annotation_custom(img) + theme(panel.background=element_blank())
}

schema_orf = function() {
    img = grid::rasterGrob(magick::image_read("external/ORF_screen_new.png"))
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
    both$pdist = with(both, sqrt(pmin(0.9,estimate.x)^2 + pmin(0.9,estimate.y)^2))

    hl = c("RBM14", "CDKN1A", "SNRPA", "ZBTB14", "POU2F1", "STAT3", "BUB1",
           "BRCA1", "RBM33", "DNMT3A", "CCNE1", "MDM2", "ERBB2", "CDC73",
           "SRSF3", "AKT3", "BAX", "POLQ", "MPP2", "ZNF879", "MCM2", "COL11A2",
           "ID2", "ZNF418", "NOMO2", "DAP3", "YY1AP1", "MSTO2P", "RFC4")
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
                   box.padding=unit(0.1, "lines"), size=3, min.segment.length=0,
                   segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        theme_classic() +
        labs(size = "TCGA\nAmplifications",
             x = "Δ Expression over expected CCLE",
             y = "Δ Expression over expected TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
}

orf_volc = function(orfdata) {
    orfdata$fill = orfdata$statistic < -5 & orfdata$estimate < log2(0.7)
    orfdata$circle =  TRUE
    plt$volcano(orfdata, label_top=35, pos_label_bias=3, max.overlaps=20) +
        labs(x = "log fold-change ORF screen",
             y = "Adjusted p-value (FDR)",
             size = "# ORFs")
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

    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
        filter(gene_name != "LOC254896") # not in tcga/ccle data

    left = (wrap_elements(schema_comp() + theme(plot.margin=margin(0,15,0,-5,"mm")))) /
        tcga_ccle_cor(comp, gistic_amp, cosmic)
    right = (wrap_elements(schema_orf()) + theme(plot.margin=margin(0,-15,0,0,"mm"))) /
        orf_volc(orfdata)

    asm = ((left + plot_layout(heights=c(1,3))) |
        (right + plot_layout(heights=c(2,3)))) + plot_layout(widths=c(3,2)) +
        plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("Fig2-Comp+ORF.pdf", 14, 10)
    print(asm)
    dev.off()
})
