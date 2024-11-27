library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

comp_orf = function(gistic_amp) {
    dset = cm$get_pancan_summary() |>
        dplyr::rename(gene_name = gene) |>
        left_join(gistic_amp) |>
        mutate(type = ifelse(is.na(type), "Background", type),
               est_ccle_tcga = (comp_ccle + comp_tcga)/2,
               label = ifelse((type == "Compensated" & stat_orf < -5) |
                    (type == "Hyperactivated" & stat_orf < -7) | stat_orf > 6 | is_argos |
                    stat_orf < -12 | abs(est_ccle_tcga) > 0.87, gene_name, NA))

    m = lm(stat_orf ~ est_ccle_tcga, data=dset) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.3f~\n~italic(P)~`=`~%.2g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)
    cols = cm$cols[c("Background", "Compensated", "Hyperactivated")]

    ggplot(dset, aes(x=est_ccle_tcga, y=stat_orf)) +
        geom_hline(yintercept=0, linetype="dashed", color="black") +
        geom_vline(xintercept=0, linetype="dashed", color="black") +
        geom_point(aes(color=type, shape=is_tox)) +
        geom_smooth(method="lm", se=FALSE) +
        scale_color_manual(values=cols, name="Compensation\nclass") +
        scale_shape_manual(values=c("TRUE"=20, "FALSE"=1), na.translate=FALSE, name="Toxic gene") +
        ggrepel::geom_label_repel(aes(label=label, color=type), size=3,
            box.padding=unit(0.1, "lines"), min.segment.length=0,
            segment.alpha=0.3, fill="#ffffffa0", label.size=NA) +
        annotate("text", y=6, x=0.6, hjust=0, label=lab, color="blue", parse=TRUE) +
        cm$theme_classic() +
        coord_cartesian(clip="off", xlim=c(-1,NA)) +
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

    p1 = plt$matrix(df, cor ~ Var1 + Var2, geom="tile") +
        scale_fill_distiller(palette="RdBu", name=ltit, limits=c(-1,1)) +
        theme(axis.title = element_blank()) +
        coord_fixed() +
        cm$text_sizes() +
        theme(axis.title = element_blank(),
              axis.text.x = element_blank()) +
        ggtitle("Pairwise")

    p2 = plt$matrix(df, cond ~ Var1 + Var2, geom="tile") +
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
        `Schukken (protein)`="Schukken_Protein", `ours`="Ours")
    ggplot(dset, aes(x=type, y=stat_orf, fill=type)) +
        geom_boxplot(outlier.shape=NA, alpha=0.9) +
        ggsignif::geom_signif(y_position=c(3.2, 4.4, 5.5, 6.6), color="black",
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0, test=t.test,
            comparisons=list(c("All genes", "Goncalves"),
                             c("All genes", "Schukken (gene)"),
                             c("All genes", "Schukken (protein)"),
                             c("All genes", "ours"))) +
        coord_cartesian(ylim=c(-7.5, 9), clip="off") +
        scale_fill_manual(values=setNames(cm$col_study, sub("\n", " ", names(cm$col_study)))) +
        labs(fill = "Study", x = "Study", y = "Δ ORF (Wald statistic)") +
        cm$theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(dset$stat_orf[dset$type=="All genes"], na.rm=TRUE),
                   linetype="dashed", color="black")
}

wgd_compare = function() {
    comp = cm$get_comp_genes(pan=TRUE)
    tox = cm$get_tox()$`Pan-Cancer` %>% filter(is_toxic) %>% pull(gene)
    argos = intersect(comp, tox)

    make_class = . %>% mutate(`Gene class` = case_when(
        gene %in% argos ~ "ARGOS",
        gene %in% comp ~ "Compensated",
        gene %in% tox ~ "Toxic",
        TRUE ~ NA_character_
    ))
    comp_comp = function(path) {
        df1 = readxl::read_xlsx(path, sheet="panWGD+")
        df2 = readxl::read_xlsx(path, sheet="panWGD-")
        both = inner_join(df1 %>% select(gene, `Compensation WGD+`=compensation),
                          df2 %>% select(gene, `Compensation WGD-`=compensation)) %>%
        mutate(label = ifelse(gene %in% c("CDKN1A", "RBM14"), gene, NA)) %>% make_class()
        m = broom::glance(lm(`Compensation WGD-` ~ `Compensation WGD+`, data=both))
        lab = sprintf("R^2~`=`~%.3f~\n~italic(P)~`=`~%.2g", m$adj.r.squared, m$p.value) %>% sub("e", "%*%10^", .)
        plt$denspt(both, aes(x=`Compensation WGD+`, y=`Compensation WGD-`), alpha=0.2) +
            geom_point(data=both[!is.na(both$`Gene class`),], aes(color=`Gene class`),
                       shape=1, alpha=0.8, stroke=1) +
            coord_cartesian(xlim=c(-1.1,1.5), ylim=c(-1.1,1.5)) +
            annotate("text", y=1.4, x=-0.8, hjust=0, label=lab, color="blue", parse=TRUE) +
#            ggrepel::geom_label_repel(aes(label=label, color=`Gene class`), segment.alpha=0.3,
#                fill="#ffffffc0", box.padding=unit(0.1, "lines"), label.size=NA) +
            guides(alpha = "none")
    }
    p1 = comp_comp("SuppData1_CCLE-comp.xlsx") + ggtitle("CCLE") + guides(fill="none")
    p2 = comp_comp("SuppData2_TCGA-comp.xlsx") + ggtitle("TCGA")

    orf1 = readRDS("../model_orf/fitsWGD.rds")
    orf = inner_join(orf1$`panWGD+` %>% select(gene=`GENE SYMBOL`, stat_wgd=statistic),
                     orf1$`panWGD-` %>% select(gene=`GENE SYMBOL`, stat_eup=statistic)) %>%
        filter(gene != "LOC254896") %>% make_class() %>%
        mutate(label = ifelse(gene %in% c("CDKN1A", "RBM14"), gene, NA))
    m = broom::glance(lm(stat_eup ~ stat_wgd, data=orf))
    lab = sprintf("R^2~`=`~%.3f~\n~italic(P)~`=`~%.2g", m$adj.r.squared, m$p.value) %>% sub("e", "%*%10^", .)
    p3 = plt$denspt(orf, aes(x=stat_eup, y=stat_wgd, alpha=0.8)) +
        geom_point(data=orf[!is.na(orf$`Gene class`),], aes(color=`Gene class`),
                   shape=1, alpha=0.8, stroke=1) +
        annotate("text", y=3.2, x=-17, hjust=0, label=lab, color="blue", parse=TRUE) +
#        ggrepel::geom_label_repel(aes(label=label, color=`Gene class`), segment.alpha=0.3,
#            fill="#ffffffc0", box.padding=unit(0.1, "lines"), label.size=NA) +
        labs(title = "ORF", x="ORF dropout WGD+", y="ORF dropout WGD-") +
        guides(alpha = "none", fill="none")

    ((p1 | p2 | p3) & cm$theme_minimal()) + plot_layout(guides="collect")
}

tissue_compare = function() {
    comp = cm$get_comp_tissue() %>%
        filter(!grepl("[pP]an", tissue)) %>%
        select(src, tissue, gene, compensation) %>%
        mutate(compensation = pmax(-1.5, pmin(2, compensation))) %>%
        tidyr::pivot_wider(names_from=src, values_from=compensation)
    m = comp %>% group_by(tissue) %>%
        summarize(res = broom::glance(lm(TCGA ~ CCLE))) %>%
        tidyr::unnest(res) %>%
        mutate(lab = sprintf("italic(P)~`=`~%.2g", p.value) %>% sub("e", "%*%10^", .))
    p1 = plt$denspt(comp, aes(x=CCLE, y=TCGA), n_tile=19, draw_pt=0) +
        scale_fill_continuous(type="viridis", trans="log1p", breaks=c(1,5,20,100,500)) +
        facet_grid(. ~ tissue) +
        geom_label(data=m, aes(label=lab), x=0, y=1.7, color="blue",
                   parse=TRUE, fill="#ffffffc0", label.size=NA)

    toxf = "../model_orf/fits_naive.xlsx"
    sheets = readxl::excel_sheets(toxf)
    tox = lapply(sheets, readxl::read_xlsx, path=toxf) %>% setNames(sheets) %>%
        bind_rows(.id="tissue") %>% filter(!grepl("pan", tissue)) %>%
        select(tissue, gene=`GENE SYMBOL`, statistic) %>%
        mutate(statistic = pmax(-10, pmin(5, statistic)))
    comp2 = comp %>% mutate(both = (CCLE+TCGA)/2) %>% inner_join(tox)
    m = comp2 %>% group_by(tissue) %>%
        summarize(res = broom::glance(lm(statistic ~ both))) %>%
        tidyr::unnest(res) %>%
        mutate(lab = sprintf("italic(P)~`=`~%.2g", p.value) %>% sub("e", "%*%10^", .))
    p2 = plt$denspt(comp2, aes(x=both, y=statistic), n_tile=19, draw_pt=0) +
        scale_fill_distiller(palette="PuOr", trans="log1p", breaks=c(1,5,20,100,500)) +
        facet_grid(. ~ tissue) +
        geom_label(data=m, aes(label=lab), x=0, y=4, color="blue",
                   parse=TRUE, fill="#ffffffc0", label.size=NA) +
        labs(x = "Mean compensation score CCLE/TCGA",
             y = "ORF (Wald st.)")

    (p1 / p2) & cm$theme_minimal()
}

sys$run({
    gistic_amp = readRDS("../data/gistic_smooth.rds")$genes %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)
    cosmic = cm$get_cosmic_annot()

    top = (wrap_plots(dens_ov()) | comp_orf(gistic_amp) | tox_implied()) +
        plot_layout(widths=c(0.75,2,1.2))
    mid = wgd_compare()
    btm = tissue_compare()

    asm = (top / mid / btm) + plot_layout(heights=c(1,1,1)) +
        plot_annotation(tag_levels='a') &
        theme(axis.text = element_text(size=10),
              legend.text = element_text(size=10),
              plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("FigS4-Overlap.pdf", 14, 15)
    print(asm)
    dev.off()
})
