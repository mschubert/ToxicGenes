library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

schema_comp = function() {
    img = grid::rasterGrob(magick::image_read("external/comp3.png"))
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
    bord = tibble(x=c(-.3,-.3), y=c(-.3,-.3), yend=c(-Inf,-.3), xend=c(-.3,-Inf))

    p = ggplot(both, aes(x=estimate.x, y=estimate.y)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_segment(data=bord, aes(x=x, y=y, xend=xend, yend=yend),
                     linetype="dotted", size=0.3, color=cm$cols[["Compensated"]]) +
        geom_segment(data=bord, aes(x=-x, y=-y, xend=-xend, yend=-yend),
                     linetype="dotted", size=0.3, color=cm$cols[["Hyperactivated"]]) +
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
                   aes(label=gene_name, color=group), max.overlaps=12,
                   box.padding=unit(0.1, "lines"), min.segment.length=0,
                   segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        cm$theme_classic() +
        labs(size = "TCGA\nAmplifications",
             x = "Compensation score CCLE",
             y = "Compensation score TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
}

orf_volc = function(orfdata) {
    orfdata$fill = orfdata$adj.p < 1e-5 & orfdata$estimate < log2(0.7)
    orfdata$circle = TRUE
    bord = tibble(x=c(-Inf,log2(0.7)), y=c(1e-5,1e-5), yend=c(1e-5,0), xend=c(log2(0.7),log2(0.7)))

    plt$volcano(orfdata, label_top=30, pos_label_bias=3, max.overlaps=15, text.size=4) +
        geom_segment(data=bord, aes(x=x, y=y, xend=xend, yend=yend),
                     linetype="dotted", size=0.3, color="grey") +
        labs(x = "log2 fold-change ORF screen",
             y = "Adjusted p-value (FDR)",
             size = "# ORFs") +
        cm$text_sizes()
}

comp_ov = function() {
    dset = cm$get_comp_tissue() %>%
        group_by(gene, src) %>%
        summarize(comp = list(c(na.omit(tissue[is_comp])))) %>%
        tidyr::unnest(comp)

    count_ov = function(ds, excl=c()) {
        ds %>% group_by(gene) %>%
        summarize(CCLE = n_distinct(comp[src=="CCLE"]),
                  TCGA = n_distinct(comp[src=="TCGA"]),
                  both = n_distinct(intersect(comp[src=="CCLE"], comp[src=="TCGA"]))) %>%
        tidyr::pivot_longer(c(CCLE, TCGA, both)) %>%
        group_by(name, value) %>%
            summarize(n = n()) %>%
        ungroup() %>%
            tidyr::complete(name, value=1:6, fill=list(n=0)) %>%
        group_by(name) %>%
            arrange(desc(value)) %>%
            mutate(n = cumsum(n)) %>%
        ungroup() %>%
        filter(value <= 6, value > 0) %>%
        mutate(name = factor(name, levels=c("TCGA", "CCLE", "both")))
    }
    pan_g = dset %>% filter(comp == "Pan-Cancer") %>% group_by(gene) %>%
        filter(all(c("CCLE", "TCGA") %in% src)) %>% pull(gene) %>% unique()
    nums_pan = dset %>% filter(comp == "Pan-Cancer") %>% count_ov()
    nums_tis = list(
        included = dset %>% filter(comp != "Pan-Cancer") %>% count_ov(pan_g),
        excluded = dset %>% filter(comp != "Pan-Cancer", ! gene %in% pan_g) %>% count_ov(pan_g)
    ) %>% bind_rows(.id="Pan-Cancer") %>% filter(`Pan-Cancer`=="included" | name == "both")

    cols = cm$cols[c("TCGA","CCLE","Compensated")]
    names(cols)[3] = "both"
    a1 = nums_pan %>% filter(name == "both", value == 1)
    a2 = nums_tis %>% filter(name == "both", value == 1)
    p1 = ggplot(nums_pan, aes(y=n, x="Pan-Cancer", fill=name)) +
        geom_col(position="dodge", alpha=0.6) +
        geom_text(data=a1, aes(label=n), x=1.3, color=cols["both"], angle=90, hjust=-0.3) +
        scale_y_continuous(trans="log1p", breaks=c(1,10,100,1000)) +
        coord_cartesian(ylim=c(0.5, NA)) +
        scale_fill_manual(values=cols) +
        cm$theme_minimal() +
        theme(axis.title.x = element_blank(),
              legend.position = "none") +
        labs(y = "Number of\ncompensated genes")
    p2 = ggplot(nums_tis, aes(x=value, y=n, color=name)) +
        geom_line(aes(group=paste(name, `Pan-Cancer`))) +
        geom_point(aes(shape=`Pan-Cancer`), size=3, alpha=0.6) +
        geom_label(data=a2, aes(label=n), color=cols["both"], hjust=0.4,
                   vjust=c(-0.3,1.6), fill="#ffffff30", label.size=NA) +
        scale_shape_manual(values=c(included=19, excluded=1)) +
        scale_y_continuous(trans="log1p", breaks=c(1,10,100,1000)) +
        scale_x_continuous(breaks=1:6) +
        coord_cartesian(ylim=c(0.5, NA)) +
        scale_color_manual(values=cols, name="Dataset") +
        cm$theme_minimal() +
        theme(axis.title.y = element_blank()) +
        labs(x="Tissue overlap") +
        plot_layout(tag_level="new")
    p1 + p2 + plot_spacer() + plot_layout(widths=c(1,5,0.5))
}

orf_ov = function(orfdata) {
    pan_g = orfdata$`Pan-Cancer` %>% filter(is_toxic) %>% pull(gene)
    tis_g = bind_rows(orfdata, .id="tissue") %>%
        filter(!grepl("Pan-Cancer", tissue), is_toxic) %>%
        pull(gene) %>% unique()

    dset = tibble(src = c("Pan-Cancer", "≥ 1 tissue"),
                  from = c(0, length(setdiff(pan_g, tis_g))),
                  to = c(length(pan_g), length(unique(c(pan_g, tis_g)))),
                  y = c(2,-2)) %>%
        mutate(src = factor(src, levels=src))
    nums = tibble(x = c(dset$from[-1], dset$to),
                  n = diff(c(dset$from, dset$to)),
                  y = c(2, 0, -2)) %>%
        mutate(x = n/2 + c(0, x[-length(x)]))

    ggplot(dset, aes(y=y, yend=y, x=from, xend=to, color=src)) +
        geom_segment(linewidth=8, alpha=0.2) +
        geom_text(data=nums, aes(x=x, label=n, y=y), color="black", inherit.aes=FALSE) +
        guides(color=guide_legend(override.aes=list(linewidth=5))) +
        scale_color_manual(values=c("coral", "steelblue")) +
        scale_x_continuous(breaks=unique(c(dset$from, dset$to))) +
        scale_y_continuous(limits=c(-15,15)) +
        labs(color = "",
             x = "Number of Toxic Genes found") +
        cm$theme_minimal() +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              plot.margin = margin(-5,0,0,0,"mm"))
}

sys$run({
    cosmic = cm$get_cosmic_annot()
    gistic = readRDS("../data/gistic_smooth.rds")$genes
    gistic_amp = gistic %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)

    ccle = readxl::read_xlsx("TableS1_CCLE-comp.xlsx", sheet="Pan-Cancer") %>%
        mutate(estimate = pmax(-2, pmin(compensation, 2.5)))
    tcga3 = readxl::read_xlsx("TableS2_TCGA-comp.xlsx", sheet="Pan-Cancer") %>%
        mutate(estimate = pmax(-2, pmin(compensation, 2.5)))
    comp_all = inner_join(ccle, tcga3, by="gene") %>%
        dplyr::rename(gene_name = gene) %>%
        left_join(cosmic)
    comp = comp_all %>% inner_join(gistic_amp)

    orfdata = sapply(readxl::excel_sheets("TableS3_ORF-toxicity.xlsx"), readxl::read_xlsx,
                     path="TableS3_ORF-toxicity.xlsx", simplify=FALSE)

    left = (wrap_elements(schema_comp() + theme(plot.margin=margin(0,0,0,-10,"mm")))) /
        tcga_ccle_cor(comp, gistic_amp, cosmic) /
        wrap_elements(comp_ov())
    right = (wrap_elements(schema_orf() + theme(plot.margin=margin(-20,-20,-15,-5,"mm")))) /
        orf_volc(orfdata$`Pan-Cancer`) /
        wrap_elements(orf_ov(orfdata))

    asm = ((left + plot_layout(heights=c(1.2,3,1.2))) |
        (right + plot_layout(heights=c(5,8,1)))) + plot_layout(widths=c(3,2)) +
        plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("Fig2-Comp+ORF.pdf", 14, 12)
    print(asm)
    dev.off()
})
