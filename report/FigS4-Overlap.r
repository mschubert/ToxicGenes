library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

comp_orf = function(all, gistic_amp) {
    dset = left_join(all %>% dplyr::rename(gene_name=gene), gistic_amp) %>%
        mutate(type = case_when(
                   est_tcga < -0.3 & est_ccle < -0.3 ~ "Compensated",
                   est_tcga > 0.3 & est_ccle > 0.3 ~ "Hyperactivated",
                   TRUE ~ "Background"
               ),
               est_ccle_tcga = (est_ccle + est_tcga)/2,
               dropout = stat_orf < -5,
               label = ifelse((type != "Background" & abs(stat_orf) > 4) | stat_orf > 6 |
                              stat_orf < -15 | abs(est_ccle_tcga) > 0.8, gene_name, NA))

    m = lm(stat_orf ~ est_ccle_tcga, data=dset) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.3f~\n~italic(P)~`=`~%.2g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    ggplot(dset, aes(x=(est_ccle+est_tcga)/2, y=stat_orf)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_point(aes(color=type, alpha=dropout)) +
        geom_smooth(method="lm", se=FALSE) +
        ggrepel::geom_label_repel(aes(label=label, color=type), size=3,
            box.padding=unit(0.1, "lines"), min.segment.length=0,
            segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        scale_color_manual(values=cm$cols[c("Background", "Compensated", "Hyperactivated")], name="Compensation\nclass") +
        scale_alpha_manual(values=c("TRUE"=0.95, "FALSE"=0.3), na.translate=FALSE, name="Toxic gene") +
        annotate("text", y=10, x=0.6, hjust=0, label=lab, color="blue", parse=TRUE) +
        cm$theme_classic() +
        coord_cartesian(clip="off") +
        labs(x = "Mean compensation score CCLE/TCGA",
             y = "ORF dropout (Wald statistic)")
}

go_cors = function() {
    orf = readxl::read_xlsx("../orf/pan/GO_Biological_Process_2021.xlsx")
    ccle_go = readRDS("../ccle/pan/stan-nb/all/GO_Biological_Process_2021.rds")$amp %>%
        select(label, stat_ccle=statistic)
    tcga_go = readRDS("../tcga/pan/stan-nb_puradj/all/GO_Biological_Process_2021.rds")[[1]] %>%
        select(label, stat_tcga=statistic, size_used)
    both = inner_join(ccle_go, tcga_go) %>% filter(size_used < 1000) %>%
        inner_join(orf %>% select(label=name, stat_orf=statistic)) %>%
        mutate(tcga_ccle = (stat_ccle+stat_tcga)/2)
    both$label[abs(both$tcga_ccle) < 1 & !grepl("33139", both$label)] = NA
    al = grep("33139", both$label, value=TRUE)

    m = broom::glance(lm(tcga_ccle ~ stat_orf, data=both))
    lab = sprintf("R^2~`=`~%.2f~\n~italic(P)~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)
    both$stat_orf = pmax(-25, both$stat_orf)

    plt$denspt(both, aes(x=tcga_ccle, y=stat_orf, label=label), size=size_used, always_label=al,
               max_ov=10, draw_label=35, palette="Greys", alpha=0.6, pal_alpha=0.5, tsize=3.5) +
        guides(alpha=FALSE) +
        scale_size_area(max_size=8, breaks=c(10,100,500,1000), name="Genes in set") +
        cm$theme_minimal() +
        annotate("text", x=6, y=-8, color="blue", label=lab, parse=TRUE) +
        labs(x = "Mean compensation score CCLE/TCGA (Wald statistic)",
             y = "Mean dropout ORF screen (Wald statistic)")
}

dens_ov = function() {
    res = readRDS("Fig3-Overlap.rds") %>%
        mutate(x = sprintf("%s:%i", chr, round(x / 1e6))) %>%
        select(-chr) %>%
        tidyr::spread(x, y)
    mat = na.omit(t(data.matrix(res[-1])))
    colnames(mat) = res$type
    ltit = "Pearson\ncorrelation\nof density"

    df = expand.grid(colnames(mat), colnames(mat)) %>%
        rowwise() %>%
        mutate(cor = cor(mat[,Var1], mat[,Var2]),
               cond = ppcor::pcor.test(mat[,Var1], mat[,Var2], mat[,"Genes"])$estimate,
               left = case_when(
                   abs(cond/cor) < 0.25 ~ "≥ 75%",
                   abs(cond/cor) < 0.6 ~ "≥ 60%"
               )
        )
    df$cond[df$Var1 == df$Var2 | df$Var1 == "Genes" | df$Var2 == "Genes"] = NA

    p1 = plt$matrix(df, cor ~ Var1 + Var2, geom="tile", text_size=10) +
        scale_fill_distiller(palette="RdBu", name=ltit, limits=c(-1,1)) +
        theme(axis.title = element_blank()) +
        coord_fixed() +
        cm$text_sizes() +
        theme(axis.title = element_blank(),
              axis.text.x = element_blank()) +
        ggtitle("Pairwise")

    p2 = plt$matrix(df, cond ~ Var1 + Var2, geom="tile", text_size=10) +
        geom_text(aes(label=left)) +
        scale_discrete_manual("label", guide=guide_legend(title="Correlation\nexplained"),
            values=c("≥ 75%"="×", "≥ 60%"="o"), na.translate=FALSE) +
        scale_fill_distiller(palette="RdBu", name=ltit, limits=c(-1,1)) +
        cm$text_sizes() +
        theme(axis.title = element_blank()) +
        coord_fixed() +
        plot_layout(tag_level="new") +
        ggtitle("Conditioned\non genes")

    (p1 / p2) + plot_layout(guides="collect")
}

tox_implied = function() {
    dset = readRDS("../misc/reviewer1/compensation.rds")$genes
    dset$type = forcats::fct_recode(dset$type, `Schukken (gene)`="Schukken_Gene",
                                    `Schukken (protein)`="Schukken_Protein")
    ggplot(dset, aes(x=type, y=stat_orf, fill=type)) +
        geom_boxplot(outlier.shape=NA, alpha=0.7) +
        ggsignif::geom_signif(y_position=c(3.2, 4.4, 5.5, 6.6), color="black", test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("All genes", "Goncalves"),
                             c("All genes", "Schukken (gene)"),
                             c("All genes", "Schukken (protein)"),
                             c("All genes", "Ours"))) +
        coord_cartesian(ylim=c(-7.5, 9), clip="off") +
        labs(fill = "Study", x = "Study", y = "Δ ORF (Wald statistic)") +
    #    scale_fill_manual(values=cm$cols[c("Background", "Compensated", "Hyperactivated")]) +
        theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(dset$stat_orf[dset$type=="All genes"], na.rm=TRUE),
                   linetype="dashed", color="black")
}

sys$run({
    gistic_amp = readRDS("../data/gistic_smooth.rds")$genes %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)
    cosmic = cm$get_cosmic_annot()
    all = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")

    top = (wrap_plots(dens_ov()) | comp_orf(all, gistic_amp)) + plot_layout(widths=c(3,4))
    mid = (tox_implied() + plot_spacer()) + plot_layout(widths=c(1,2))

    asm = (top / mid) + plot_layout(heights=c(3,2)) +
        plot_annotation(tag_levels='a') &
        theme(axis.text = element_text(size=10),
              legend.text = element_text(size=10),
              plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("FigS4-Overlap.pdf", 11, 9)
    print(asm)
    dev.off()
})
