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
                              stat_orf < -12 | abs(est_ccle_tcga) > 0.8, gene_name, NA))

    m = lm(stat_orf ~ est_ccle_tcga, data=dset) %>% broom::glance()
    lab = sprintf("R^2~`=`~%.3f~\n~p~`=`~%.2g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    ggplot(dset, aes(x=(est_ccle+est_tcga)/2, y=stat_orf)) +
        geom_hline(yintercept=0, size=2, linetype="dashed", color="grey") +
        geom_vline(xintercept=0, size=2, linetype="dashed", color="grey") +
        geom_point(aes(color=type, alpha=dropout)) +
        ggrepel::geom_label_repel(aes(label=label, color=type), size=3,
            box.padding=unit(0.1, "lines"), min.segment.length=0,
            segment.alpha=0.3, fill="#ffffff50", label.size=NA) +
        scale_color_manual(values=cm$cols[c("Background", "Compensated", "Hyperactivated")]) +
        scale_alpha_manual(values=c("TRUE"=0.95, "FALSE"=0.4), name="Dropout") +
        annotate("text", y=10, x=0.6, hjust=0, label=lab, color="blue", parse=TRUE) +
        geom_smooth(method="lm", se=FALSE) +
        theme_classic() +
        coord_cartesian(clip="off") +
        labs(x = "Mean expression over expected CCLE/TCGA",
             y = "ORF log fold-chance (Wald statistic)")
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

    m = broom::glance(lm(tcga_ccle ~ stat_orf, data=both))
    lab = sprintf("R^2~`=`~%.2f~\n~p~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    plt$denspt(both, aes(x=tcga_ccle, y=stat_orf, label=label), size=size_used,
               max_ov=20, draw_label=120) +
        scale_size_area(max_size=8, breaks=c(10,100,500,1000), name="Genes in set") +
        theme_minimal() +
        ylim(c(-24,NA)) +
        annotate("text", x=6, y=-11, color="blue", label=lab, parse=TRUE) +
        labs(x = "Mean expression over expected CCLE/TCGA (Wald statistic)",
             y = "Mean dropout ORF screen (Wald statistic)")
}

sys$run({
    gistic_amp = readRDS("../data/gistic_smooth.rds")$genes %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene_name, frac)
    cosmic = cm$get_cosmic_annot()
    all = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")

    top = (comp_orf(all, gistic_amp) | plot_spacer()) + plot_layout(widths=c(3,2))
    mid = go_cors()

    asm = (top / mid) + plot_layout(heights=c(1,2)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS4-Overlap.pdf", 11, 14)
    print(asm)
    dev.off()
})
