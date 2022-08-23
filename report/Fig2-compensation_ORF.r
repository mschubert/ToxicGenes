library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
fig1 = import('./Fig1-Motivation')

tcga_ccle_cor = function(gistic_amp, cosmic) {
    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    both = inner_join(ccle, tcga3, by="gene") %>%
        dplyr::rename(gene_name = gene) %>%
        left_join(cosmic) %>%
        inner_join(gistic_amp)

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

    m = lm(estimate.y ~ estimate.x, data=both) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.2f~\n~p~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    p = ggplot(both, aes(x=estimate.x, y=estimate.y)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_rect(xmin=-Inf, ymin=-Inf, xmax=-0.3, ymax=-0.3, color="orange",
                  linetype="dashed", aes(fill="Compensated")) +
        scale_fill_manual(values=c(Compensated="#fff9f5de"), name="Area") +
        geom_point(data=both %>% filter(is.na(type)),
                   aes(size=n_aneup.y), alpha=0.2) +
        geom_density2d(color="white", breaks=c(0.5, 0.99)) +
        geom_point(data=both %>% filter(!is.na(type)),
                   aes(size=n_aneup.y, color=type), alpha=0.7) +
        geom_smooth(method="lm", color="blue") +
        annotate("text", y=1.3, x=-0.9, hjust=0, label=lab, color="blue", parse=TRUE) +
        ggrepel::geom_label_repel(data=both %>% filter(!is.na(type)),
                   aes(label=gene_name, color=type), max.overlaps=11,
                   size=3, min.segment.length=0,
                   segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        theme_classic() +
        labs(x = "Expression over expected CCLE",
             y = "Expression over expected TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
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
    cosmic = fig1$get_cosmic_annot()
    both = left_join(orfdata, cosmic) %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = ifelse(type == "Both", "OG+TSG", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG", "OG+TSG")))

    ggplot(both, aes(x=type, y=statistic, fill=type)) +
        geom_boxplot(outlier.shape=NA) +
        ggsignif::geom_signif(y_position=c(6.5, 9), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))) +
        coord_cartesian(ylim=c(-8, 11)) +
        labs(fill = "Driver status", x = "Gene type subset", y = "Δ ORF (Wald statistic)") +
        theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

amp_del_orf = function(orfdata) {
    gwide = fig1$get_gistic_scores() %>%
        tidyr::pivot_wider(names_from="type", values_from="frac") %>%
        mutate(type = case_when(
            amplification > 0.15 & deletion < -0.15 ~ "Amp+Del",
            amplification > 0.15 ~ "Amplified",
            deletion < -0.15 ~ "Deleted",
            TRUE ~ "Background"
        ))

    both = inner_join(orfdata, gwide) %>%
        mutate(type = factor(type, levels=c("Background", "Amplified", "Deleted", "Amp+Del")))

    ggplot(both, aes(x=type, y=statistic, fill=type)) +
        geom_boxplot(outlier.shape=NA) +
        ggsignif::geom_signif(y_position=c(5, 6.5), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))) +
        coord_cartesian(ylim=c(-5, 8)) +
        labs(fill = "Frequent CNA", x = "Copy number subset", y = "Δ ORF (Wald statistic)") +
        theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

sys$run({
    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
        filter(gene_name != "LOC254896") # not in tcga/ccle data

    gistic_amp = fig1$get_gistic_scores() %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)
    cosmic = fig1$get_cosmic_annot()

    orf_cors = (amp_del_orf(orfdata) | og_tsg_orf(orfdata))

    top = tcga_ccle_cor(gistic_amp, cosmic) | (plt$text("comp GO") / plt$text("TCGA/CCLE OG/TSG"))
    btm = (orf_volc(orfdata) | (plt$text("orf GO") / orf_cors)) + plot_layout(widths=c(1,1.2))

    asm = (top / btm) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("Fig2-compensation_ORF.pdf", 16, 14)
    print(asm)
    dev.off()
})
