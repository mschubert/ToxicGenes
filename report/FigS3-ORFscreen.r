library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
orf = import('../orf/overview_naive')
cm = import('./common')

facet_plot = function(ov, mapping, hl=c("RBM14", "CDKN1A")) {
#    loess_sd = ov %>% group_by(cells) %>%
#        summarize(mod = list(msir::loess.sd(!! aes$x, !! aes$y))) %>%
#        rowwise() %>%
#        mutate(res = list(data.frame(x=mod$x, y=mod$y, sd=mod$sd))) %>%
#        select(-mod) %>%
#        tidyr::unnest(res)

    hl = ov %>% dplyr::rename(Gene=`GENE SYMBOL`) %>% filter(Gene %in% hl)

    ggplot(ov, mapping) +
        geom_hex(aes(color=..count.., fill=..count..), bins=50) +
        geom_hline(yintercept=0, color="grey", linetype="dashed") +
        scale_fill_continuous(type = "viridis", trans="log1p", breaks=c(1,5,20,100,500)) +
        scale_color_continuous(type = "viridis", trans="log1p", guide="none") +
        ggnewscale::new_scale(c("color", "fill")) +
        geom_point(data=hl, aes(color=Gene)) +
#        geom_line(data=loess_sd, aes(x=x, y=y), color="red") +
#        geom_line(data=loess_sd, aes(x=x, y=y+sd), color="red", linetype="dashed") +
#        geom_line(data=loess_sd, aes(x=x, y=y-sd), color="red", linetype="dashed") +
        facet_wrap(~ cells, ncol=6) +
        theme_minimal() +
        theme(strip.background = element_rect(color=NA, fill="#ffffffc0"))
}

screen_cor = function(ov) {
    wide = ov %>%
        select(`Construct IDs`, cells, `LFC DMSO/ETP`) %>%
        tidyr::pivot_wider(names_from=cells, values_from=`LFC DMSO/ETP`)
    mat = data.matrix(wide[-1])
    rownames(mat) = wide[[1]]

    cmat = cor(mat) %>%
        reshape2::melt() %>%
        plt$cluster(value ~ Var1 + Var2)
    plt$matrix(cmat, value ~ Var1 + Var2, geom="tile") +
        scale_fill_distiller(palette="RdBu", name="Pearson\ncorrelation") +
        theme(axis.title = element_blank()) +
        coord_fixed() +
        theme(axis.text.x = element_text(angle=90, hjust=1))
}

go_volc = function() {
    res = readxl::read_xlsx("../orf/pan/GO_Biological_Process_2018.xlsx")
    plt$volcano(res, label_top=35) + guides(size="none") +
        xlab("Mean z-score LFC")
}

og_tsg_orf = function(gistic, orfdata) {
    freq_amp_genes = gistic %>%
        filter(type == "amplification", frac > 0.15) %>% pull(gene_name)
    cosmic = cm$get_cosmic_annot()
    both = left_join(orfdata, cosmic) %>%
        filter(gene_name %in% freq_amp_genes) %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG"))) %>%
        filter(!is.na(type))

    ggplot(both, aes(x=type, y=statistic, fill=type)) +
        geom_boxplot(outlier.shape=NA, alpha=0.7) +
        ggsignif::geom_signif(y_position=c(6.5, 9), color="black", test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))) +
        coord_cartesian(ylim=c(-8, 11)) +
        scale_fill_manual(values=cm$cols[c("Background", "Oncogene", "TSG")]) +
        labs(fill = "Driver status\n(freq. amplified)", x = "Gene type subset", y = "Δ ORF (Wald statistic)") +
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
        mutate(type = factor(type, levels=c("Background", "Amplified", "Deleted"))) %>%
        filter(!is.na(type))

    ggplot(both, aes(x=type, y=statistic, fill=type)) +
        geom_boxplot(outlier.shape=NA, alpha=0.7) +
        ggsignif::geom_signif(y_position=c(5, 6.5), color="black", test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))) +
        coord_cartesian(ylim=c(-5, 8)) +
        scale_fill_manual(values=cm$cols[c("Background", "Amplified", "Deleted")]) +
        labs(fill = "Frequent CNA", x = "Copy number subset", y = "Δ ORF (Wald statistic)") +
        theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")$genes
    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
        filter(gene_name != "LOC254896") # not in tcga/ccle data
    ov = readRDS("../orf/overview.rds") %>%
        mutate(cells = sprintf("%s (%s)", cells, tissue))

    naive = facet_plot(ov, aes(x=DMSO, y=`LFC DMSO/ETP`)) +
        ylim(c(-4,4)) +
        coord_cartesian(ylim=c(-2.6,2.6), clip="off")

    left = (naive / (og_tsg_orf(gistic, orfdata) | amp_del_orf(gistic, orfdata))) +
        plot_layout(heights=c(2,1))
    right = (wrap_elements(screen_cor(ov)) / go_volc()) +
        plot_layout(heights=c(1,2))

    asm = (left | right) + plot_layout(widths=c(3,2)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("FigS3-ORFscreen.pdf", 14, 10)
    print(asm)
    dev.off()
})
