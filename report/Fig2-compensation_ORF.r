library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

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

    comphyp = data.frame(
        xmin = c(-Inf, 0.3), xmax = c(-0.3, Inf),
        ymin = c(-Inf, 0.3), ymax = c(-0.3, Inf),
        fill = c("Compensated", "Hyperactivated")
    )

    m = lm(estimate.y ~ estimate.x, data=both) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.2f~\n~p~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    p = ggplot(both, aes(x=estimate.x, y=estimate.y)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_rect(data=comphyp, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            color=cm$cols[comphyp$fill], fill="#ffffff00", linetype="dashed", inherit.aes=FALSE) +
        geom_point(data=both %>% filter(is.na(type)),
                   aes(size=n_aneup.y), alpha=0.2) +
        geom_density2d(color="white", breaks=c(0.5, 0.99)) +
        geom_point(data=both %>% filter(!is.na(type)),
                   aes(size=n_aneup.y, color=type), alpha=0.7) +
        scale_color_manual(values=cm$cols[c("Oncogene", "TSG", "OG+TSG")], name="Driver status") +
        geom_smooth(method="lm", color="blue") +
        annotate("text", x=0.2, y=-0.8, label="Compensated", color=cm$cols[["Compensated"]],
                 size=6, fontface="bold", alpha=0.7) +
        annotate("text", x=0.4, y=1.55, label="Hyperactivated", color=cm$cols[["Hyperactivated"]],
                 size=4, fontface="bold", alpha=0.7, hjust=0) +
        annotate("text", y=1.3, x=-0.9, hjust=0, label=lab, color="blue", parse=TRUE) +
        ggrepel::geom_label_repel(data=both %>% filter(!is.na(type)),
                   aes(label=gene_name, color=type), max.overlaps=11,
                   size=3, min.segment.length=0,
                   segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        theme_classic() +
        labs(size = "Aneuploid\nsamples",
             x = "Expression over expected CCLE",
             y = "Expression over expected TCGA") +
        plot_layout(tag_level="new")

    dx + plot_spacer() + plot_spacer() + p + dy + guide_area() +
        plot_layout(widths=c(10,1,2), heights=c(1,10), guides="collect")
}

go_tcga_ccle = function() {
    ccle = readRDS("../ccle/pan/stan-nb/all/GO_Biological_Process_2021.rds")$amp %>%
        transmute(label=label, est_CCLE=estimate, fdr_CCLE=adj.p)
    tcga = readRDS("../tcga/pan/stan-nb_puradj/all/GO_Biological_Process_2021.rds")[[1]] %>%
        transmute(label=label, est_TCGA=estimate, fdr_TCGA=adj.p)
    dset = inner_join(ccle, tcga) %>%
        mutate(mean_est = rowMeans(cbind(est_CCLE, est_TCGA))) %>%
        arrange(mean_est) %>%
        head(15) %>%
        mutate(label = forcats::fct_reorder(label, mean_est, .desc=TRUE)) %>%
        tidyr::pivot_longer(c(est_CCLE, est_TCGA),
                            names_pattern="_(TCGA|CCLE)",
                            names_to="dset")
    txt = dset %>% select(label) %>% distinct()

    ggplot(dset, aes(x=label, y=value)) +
        geom_col(aes(fill=dset), width=1.2, position=position_dodge(width=0.6), alpha=0.5) +
        geom_text(data=txt, aes(label=paste(" ", label)), y=0, hjust=0) +
        coord_flip(expand=FALSE, clip="off") +
        scale_y_reverse() +
        scale_fill_manual(values=cm$cols[c("TCGA", "CCLE")], name="Dataset") +
        labs(x = "Gene Ontology", y = "Mean compensation score in CCLE/TCGA") +
        theme_classic() +
        theme(axis.text.y = element_blank())
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

comp_tcga_ccle = function(comp) {
    both = comp %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG", "OG+TSG")))

    common = function(y, coordy, sigy) list(
        geom_boxplot(outlier.shape=NA, alpha=0.7),
        ggsignif::geom_signif(y_position=sigy, color="black", test=wilcox.test,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))),
        scale_fill_manual(values=cm$cols[c("Background", "Oncogene", "TSG", "OG+TSG")]),
        labs(fill = "Driver status", x = "Gene type subset"),
        theme_classic(),
        coord_cartesian(ylim=coordy),
        theme(axis.text.x = element_blank()),
        geom_hline(yintercept=median(y[both$type=="Background"]),
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
        ggsignif::geom_signif(y_position=c(6.5, 9), color="black", test=wilcox.test,
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
        geom_boxplot(outlier.shape=NA, alpha=0.7) +
        ggsignif::geom_signif(y_position=c(5, 6.5), color="black", test=wilcox.test,
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
    gistic_amp = gistic %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)

    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    comp = inner_join(ccle, tcga3, by="gene") %>%
        dplyr::rename(gene_name = gene) %>%
        left_join(cosmic) %>%
        inner_join(gistic_amp)

    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
        filter(gene_name != "LOC254896") # not in tcga/ccle data
    orf_cors = amp_del_orf(gistic, orfdata) | og_tsg_orf(orfdata)

    top = (tcga_ccle_cor(comp, gistic_amp, cosmic) |
        ((comp_tcga_ccle(comp) / go_tcga_ccle()) + plot_layout(heights=c(2,3)))) +
        plot_layout(widths=c(1.3,1))
    btm = (orf_volc(orfdata) | ((orf_cors / go_orf()) + plot_layout(heights=c(2,3)))) +
        plot_layout(widths=c(1,1.3))

    asm = (top / btm) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("Fig2-compensation_ORF.pdf", 15, 14)
    print(asm)
    dev.off()
})