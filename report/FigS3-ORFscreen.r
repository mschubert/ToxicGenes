library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
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
        geom_hex(aes(fill=..count..), bins=30) +
        geom_hline(yintercept=0, color="grey", linetype="dashed") +
        scale_fill_continuous(type="viridis", trans="log1p", breaks=c(1,5,20,100,500)) +
        ggnewscale::new_scale(c("fill")) +
        geom_point(data=hl, aes(fill=Gene), color="black", shape=21, size=2) +
#        geom_line(data=loess_sd, aes(x=x, y=y), color="red") +
#        geom_line(data=loess_sd, aes(x=x, y=y+sd), color="red", linetype="dashed") +
#        geom_line(data=loess_sd, aes(x=x, y=y-sd), color="red", linetype="dashed") +
        facet_wrap(~ cells, ncol=6) +
        cm$theme_minimal() +
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
    levels(cmat$Var1) = levels(cmat$Var2)
    plt$matrix(cmat, value ~ Var1 + Var2, geom="tile") +
        scale_fill_distiller(palette="RdBu", name="Pearson\ncorrelation") +
        coord_fixed() +
        cm$text_sizes() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
}

go_volc = function() {
    res = readxl::read_xlsx("../model_orf/pan/GO_Biological_Process_2018.xlsx")
    res$name = sub("positive", "pos.", res$name)
    res$name = sub("regulation", "reg.", res$name)
    res$name = sub("activation", "act.", res$name)
    res$name = sub("phosphorylation", "phosph.", res$name)
    res$name = sub("natural killer", "NK", res$name)
    res$name = sub("-mediated", "", res$name)
    res$name = sub("peptidyl-serine", "pept-ser.", res$name)
    res$name = sub("involved", "inv.", res$name)
    res$name = sub("response", "resp.", res$name)
    res$name[abs(res$estimate) < 0.5 | res$adj.p > 1e-30] = NA
    res$name[grepl("GO:1990000|GO:00(42100|35025|30098|30857|34393|33139|43331)", res$name)] = NA # labels overlap otherwise
    plt$volcano(res, label_top=35, max.overlaps=10, text.size=3.5) + guides(size="none") +
        xlab("Mean z-score LFC") +
        ylab("Upper bound FDR (values too close to zero)") +
        coord_cartesian(ylim=c(1,1e-95)) +
        cm$text_sizes()
}

og_tsg_orf = function(gistic, orfdata) {
    freq_amp_genes = gistic %>%
        filter(type == "amplification", frac > 0.15) %>% pull(gene)
    cosmic = cm$get_cosmic_annot() %>% dplyr::rename(gene=gene_name)
    both = left_join(orfdata, cosmic) %>%
        filter(gene %in% freq_amp_genes) %>%
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
        labs(fill = "Driver status\n(freq. amplified)", x = "Driver status", y = "Δ ORF (Wald statistic)") +
        cm$theme_classic() +
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
        ggsignif::geom_signif(y_position=c(4.7, 6.5), color="black", test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))) +
        coord_cartesian(ylim=c(-5, 8.5)) +
        scale_fill_manual(values=cm$cols[c("Background", "Amplified", "Deleted")]) +
        labs(fill = "Frequent CNA", x = "Frequent CNA", y = "Δ ORF (Wald statistic)") +
        cm$theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(both$statistic[both$type=="Background"]),
                   linetype="dashed", color="black")
}

tissue_ov = function(orfdata) {
    meta = readr::read_tsv("../data/orf/tissues.txt")
    dset = bind_rows(orfdata, .id="screen") %>%
        left_join(meta %>% select(screen=cells, tissue)) %>%
        mutate(tissue = ifelse(screen=="Pan-Cancer", "Pan-Cancer", tissue),
               tissue = factor(tissue) %>% relevel("Pan-Cancer"),
               screen = ifelse(screen=="Pan-Cancer", "all", screen)) %>%
        group_by(gene) %>%
            filter(sum(is_toxic) >= 6) %>%
        ungroup() %>%
        mutate(s = ifelse(is_toxic, 1, 0.7),
               estimate = pmax(-3, pmin(3, estimate)))

    ggplot(dset, aes(x=gene, y=forcats::fct_rev(screen), fill=estimate)) +
        geom_tile(aes(width=s, height=s)) +
        scale_fill_distiller(palette="PuOr", limits=max(abs(dset$estimate))*c(-1,1), name="log2 FC") +
        scale_y_discrete(position = "right") +
        facet_grid(tissue ~ ., scales="free", space="free", switch="y") +
        cm$theme_minimal() +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),
              strip.text.y.left = element_text(angle=0, hjust=0),
              strip.background = element_rect(color="black", linewidth=1)) +
        labs(x = "Gene",
             y = "Cell line")
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")$genes %>% dplyr::rename(gene=gene_name)
    orfdata = cm$get_tox()
    ov = readRDS("../model_orf/overview.rds")

    naive = facet_plot(ov, aes(x=DMSO, y=`LFC DMSO/ETP`)) +
        ylim(c(-4,4)) +
        xlab("Log10 read count DMSO condition") +
        coord_cartesian(ylim=c(-2.6,2.6), clip="off")

    left = (naive / (og_tsg_orf(gistic, orfdata$`Pan-Cancer`) |
                     amp_del_orf(gistic, orfdata$`Pan-Cancer`))) +
        plot_layout(heights=c(2,1))
    right = (wrap_elements(go_volc() + theme(plot.margin=margin(-10,0,0,0,"mm"))) /
             wrap_elements(screen_cor(ov) + theme(plot.margin=margin(-5,0,-15,-5,"mm")))) +
        plot_layout(heights=c(2,1.2))
    btm = tissue_ov(orfdata)

    asm = (((left | right) + plot_layout(widths=c(3,2))) / wrap_elements(btm)) +
        plot_layout(heights=c(10,6)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("FigS3-ORFscreen.pdf", 14, 14)
    print(asm)
    dev.off()
})
