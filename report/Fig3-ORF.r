library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

schema = function() {
    img = grid::rasterGrob(magick::image_read("external/tox.svg"))
    ggplot() + annotation_custom(img) + theme(panel.background=element_blank())
}

go_orf = function() {
    dset = readxl::read_xlsx("../orf/pan/GO_Biological_Process_2018.xlsx") %>%
        filter(adj.p < 0.05) %>%
        arrange(estimate) %>%
        head(15) %>%
        mutate(label = forcats::fct_reorder(name, estimate, .desc=TRUE))

    ggplot(dset, aes(x=label, y=estimate)) +
        geom_col(aes(fill="ORF"), alpha=0.5) +
        scale_fill_manual(values=cm$cols["ORF"], name="Dataset") +
        geom_text(aes(label=paste(" ", label)), y=0, hjust=0) +
        coord_flip(expand=FALSE, clip="off") +
        scale_y_reverse() +
        labs(x = "Gene Ontology", y = "ORF dropout (mean Wald statistic)") +
        theme_classic() +
        theme(axis.text.y = element_blank())
}

orf_volc = function(orfdata) {
    orfdata$fill = orfdata$statistic < -5
    orfdata$circle =  TRUE
    plt$volcano(orfdata, label_top=35, pos_label_bias=3, max.overlaps=20) +
        labs(x = "log fold-change ORF screen",
             y = "Adjusted p-value (FDR)",
             size = "# ORFs")
}

og_tsg_orf = function(orfdata) {
    cosmic = cm$get_cosmic_annot()
    both = left_join(orfdata, cosmic) %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = ifelse(type == "Both", "OG+TSG", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG", "OG+TSG")))

    ggplot(both, aes(x=type, y=statistic, fill=type)) +
        geom_boxplot(outlier.shape=NA, alpha=0.7) +
        ggsignif::geom_signif(y_position=c(6.5, 9), color="black", test=t.test,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))) +
        coord_cartesian(ylim=c(-8, 11)) +
        scale_fill_manual(values=cm$cols[c("Background", "Oncogene", "TSG", "OG+TSG")]) +
        labs(fill = "Driver status", x = "Gene type subset", y = "Δ ORF (Wald statistic)") +
        theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

amp_del_orf = function(gistic, orfdata) {
    gwide = gistic %>%
        tidyr::pivot_wider(names_from="type", values_from="frac")
    gwide = gwide %>%
        mutate(type = case_when(
            amplification > 0.15 & deletion < -0.15 ~ "Amp+Del",
            amplification > 0.15 ~ "Amplified",
            deletion < -0.15 ~ "Deleted",
            TRUE ~ NA_character_
        )) %>% filter(!is.na(type)) %>% bind_rows(gwide %>% mutate(type="Background"))

    both = inner_join(orfdata, gwide) %>%
        mutate(type = factor(type, levels=c("Background", "Amplified", "Deleted", "Amp+Del")))

    ggplot(both, aes(x=type, y=statistic, fill=type)) +
        geom_boxplot(outlier.shape=NA, alpha=0.7) +
        ggsignif::geom_signif(y_position=c(5, 6.5), color="black", test=t.test,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))) +
        coord_cartesian(ylim=c(-5, 8)) +
        scale_fill_manual(values=cm$cols[c("Background", "Amplified", "Deleted", "Amp+Del")]) +
        labs(fill = "Frequent CNA", x = "Copy number subset", y = "Δ ORF (Wald statistic)") +
        theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

sys$run({
    cosmic = cm$get_cosmic_annot()
    gistic = readRDS("../data/gistic_smooth.rds")$genes

    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
        filter(gene_name != "LOC254896") # not in tcga/ccle data
    orf_cors = amp_del_orf(gistic, orfdata) | og_tsg_orf(orfdata)

    top = schema()
    btm = (orf_volc(orfdata) | ((orf_cors / go_orf()) + plot_layout(heights=c(2,3)))) +
        plot_layout(widths=c(1,1.3))

    asm = (top / btm) + plot_layout(heights=c(1,3)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("Fig3-ORF.pdf", 15, 10)
    print(asm)
    dev.off()
})
