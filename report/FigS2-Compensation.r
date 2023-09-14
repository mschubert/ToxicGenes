library(dplyr)
library(ggplot2)
library(patchwork)
`%$%` = magrittr::`%$%`
sys = import('sys')
seq = import('seq')
plt = import('plot')
gset = import('genesets')
cm = import('./common')

tcga_vs_ccle = function(hl=c("RBM14", "CDKN1A")) {
    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga1 = readxl::read_xlsx("../tcga/pan/stan-nb_naive.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga2 = readxl::read_xlsx("../tcga/pan/stan-nb_pur.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))

    dset = ccle %>% select(gene, CCLE=estimate) %>%
        left_join(tcga1 %>% select(gene, `No purity correction`=estimate)) %>%
        left_join(tcga2 %>% select(gene, `Common purity term`=estimate)) %>%
        left_join(tcga3 %>% select(gene, `Purity per tissue`=estimate)) %>%
        tidyr::gather("type", "value", -gene, -CCLE) %>%
        mutate(type = factor(type) %>% relevel("No purity correction"))
    hl = dset %>% dplyr::rename(Gene=gene) %>% filter(Gene %in% hl)

    mods = dset %>% group_by(type) %>%
        summarize(mod = list(lm(value ~ CCLE))) %>%
        rowwise() %>%
        mutate(tidy = list(broom::tidy(mod)),
               glance = list(broom::glance(mod)),
               intcp = tidy$estimate[tidy$term == "(Intercept)"],
               slope = tidy$estimate[tidy$term == "CCLE"],
               angle = atan(slope) * 180/pi) %>%
        select(-tidy, -mod) %>%
        tidyr::unnest(glance) %>%
        mutate(label = sprintf("R^2~`=`~%.2f~italic(P)~`=`~10^%.0f", adj.r.squared, ceiling(log10(p.value))))

    ggplot(dset, aes(x=CCLE, y=value)) +
        geom_vline(xintercept=0, color="grey", linetype="dashed", size=1) +
        geom_hline(yintercept=0, color="grey", linetype="dashed", size=1) +
        geom_hex(bins=50) +
        scale_fill_continuous(type = "viridis", trans="log1p", breaks=c(1,5,20,100,500)) +
        facet_wrap(~ type) +
        geom_smooth(method="lm", color="blue", se=FALSE, size=0.7) +
        geom_label(data=mods, aes(x=1, y=-1, label=label), parse=TRUE, color="blue",
                   fill="#ffffffc0", label.size=NA, hjust=1, vjust=0.8, size=4) +
        ggnewscale::new_scale(c("fill")) +
        geom_point(data=hl, aes(fill=Gene), color="black", shape=21, size=2.5) +
        labs(y = "Expression over expected TCGA",
             x = "Expression over expected CCLE") +
        coord_cartesian(ylim=c(-1.2, 1.2)) +
        theme_minimal() +
        theme(strip.text = element_text(size=12),
              strip.background = element_rect(color=NA, fill="#ffffffc0"))
}

go_cors = function() {
    ccle_go = readRDS("../ccle/pan/stan-nb/all/GO_Biological_Process_2021.rds")$amp %>%
        select(label, stat_ccle=statistic)
    tcga_go = readRDS("../tcga/pan/stan-nb_puradj/all/GO_Biological_Process_2021.rds")[[1]] %>%
        select(label, stat_tcga=statistic, size_used)
    both = inner_join(ccle_go, tcga_go) %>% filter(size_used < 1000)

    m = broom::glance(lm(stat_ccle ~ stat_tcga, data=both))
    lab = sprintf("R^2~`=`~%.2f~\n~italic(P)~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    plt$denspt(both, aes(x=stat_tcga, y=stat_ccle, label=label), size=size_used,
               palette="Greys", alpha=0.6, pal_alpha=0.5, tsize=3.5, max_ov=3) +
        scale_size_area(max_size=8, breaks=c(10,100,500,1000), name="Genes in set") +
        theme_minimal() +
        guides(alpha="none") +
        labs(title = "Gene Ontology: Biological Process",
             x = "Δ Expression over expected TCGA (Wald stat.)",
             y = "Δ Expression over expected CCLE (Wald stat.)") +
        annotate("text", x=11, y=-6.5, color="blue", label=lab, size=4, parse=TRUE)
}

cna_comp = function(gistic, comp_all) {
    gwide = gistic %>%
        tidyr::pivot_wider(names_from="type", values_from="frac")
    gwide = gwide %>%
        mutate(type = case_when(
            amplification > 0.15 & deletion < -0.15 ~ "Amp+Del",
            amplification > 0.15 ~ "Amplified",
            deletion < -0.15 ~ "Deleted",
            TRUE ~ NA_character_
        )) %>% filter(!is.na(type)) %>% bind_rows(gwide %>% mutate(type="Background"))

    both = inner_join(comp_all %>% select(-type), gwide) %>%
        mutate(type = factor(type, levels=c("Background", "Amplified", "Deleted"))) %>%
        filter(!is.na(type))

    common = function(y, coordy, sigy) list(
        geom_boxplot(outlier.shape=NA, alpha=0.7),
        ggsignif::geom_signif(y_position=sigy, color="black", test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Background", "Amplified"), c("Background", "Deleted"))),
        scale_fill_manual(values=cm$cols[c("Background", "Amplified", "Deleted")]),
        labs(fill = "Frequent CNA", x = "Copy number subset", y = "Δ ORF (Wald statistic)"),
        theme_classic(),
        coord_cartesian(ylim=coordy, clip="off"),
        theme(axis.text.x = element_blank()),
        geom_hline(yintercept=median(y[both$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
    )

    p1 = ggplot(both, aes(x=type, y=estimate.x, fill=type)) +
        common(both$estimate.x, c(-0.3, 0.7), c(0.4, 0.55)) +
        labs(title = "CCLE", y="Δ Expression / expected")
    p2 = ggplot(both, aes(x=type, y=estimate.y, fill=type)) +
        common(both$estimate.y, c(-0.5, 1.4), c(0.9, 1.2)) +
        labs(title = "TCGA", y="")

    (p1 | (p2 + plot_layout(tag_level="new"))) + plot_layout(guides="collect")
}

og_comp = function(comp) {
    both = comp %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG"))) %>%
        filter(!is.na(type))

    common = function(y, coordy, sigy) list(
        geom_boxplot(outlier.shape=NA, alpha=0.7),
        ggsignif::geom_signif(y_position=sigy, color="black", test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Background", "Oncogene"), c("Background", "TSG"))),
        scale_fill_manual(values=cm$cols[c("Background", "Oncogene", "TSG")]),
        labs(fill = "Driver status\n(freq. amplified)", x = "Gene type subset"),
        theme_classic(),
        coord_cartesian(ylim=coordy, clip="off"),
        theme(axis.text.x = element_blank()),
        geom_hline(yintercept=median(y[both$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
    )

    p1 = ggplot(both, aes(x=type, y=estimate.x, fill=type)) +
        common(both$estimate.x, c(-0.3, 0.7), c(0.38, 0.55)) +
        labs(title = "CCLE", y="Δ Expression / expected")
    p2 = ggplot(both, aes(x=type, y=estimate.y, fill=type)) +
        common(both$estimate.y, c(-0.5, 1.4), c(0.9, 1.2)) +
        labs(title = "TCGA", y="")

    (p1 | (p2 + plot_layout(tag_level="new"))) + plot_layout(guides="collect")
}

rpe_comp = function(all) {
    lookup = c(SS6="chr +7", SS51="+7 +22", SS111="+8 +9 +18")
    chrs = seq$gene_table() %>% select(label=external_gene_name, chr=chromosome_name) %>% distinct()
    comp = all %>% filter(est_ccle < -0.3, est_tcga < -0.3) %>% pull(gene)
    dset = readRDS("../data/rnaseq_rpe1_broad/compute_fcs.rds") %>%
        tidyr::unnest(genes) %>%
        inner_join(chrs) %>%
        filter((term == "SS6" & chr == "7") |
               (term == "SS51" & chr %in% c("7", "22")) |
               (term == "SS111" & chr %in% c("8", "9", "18"))) %>%
        transmute(Sample = factor(lookup[term], levels=lookup),
                  Gene=label, chr=chr, LFC=log2FoldChange,
                  status = ifelse(Gene %in% comp, "Compensated", "Background"),
                  status = factor(status, levels=c("Background", "Compensated"))) %>%
        group_by(Sample, chr) %>%
            mutate(LFC = scale(LFC, scale=FALSE)[,1]) %>%
        ungroup()

    ggplot(dset, aes(x=status, y=2^LFC, color=status)) +
        geom_boxplot(aes(fill=status), outlier.shape=NA, alpha=0.3) +
        ggbeeswarm::geom_quasirandom(dodge.width=0.8, aes(alpha=status)) +
        scale_y_log10() +
        facet_wrap(~ Sample) +
        coord_cartesian(ylim=c(0.4, 2)) +
        labs(title = "Isogenic RPE-1 lines",
             x = "Clone with chromosome amplification",
             y = "Fold-change amplified chr\nvs. whole chromosomes") +
        scale_color_manual(values=c(cm$cols[c("Background", "Compensated")]), name="Genes") +
        scale_fill_manual(values=c(cm$cols[c("Background", "Compensated")]), name="Genes") +
        scale_alpha_manual(values=c(Background=0.1, Compensated=0.6), guide="none") +
        theme_minimal() +
        theme(axis.text.x = element_blank(),
              strip.text = element_text(size=12)) +
        ggsignif::geom_signif(color="black", y_position=-0.15,
            test=function(...) t.test(..., alternative="greater"),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Background", "Compensated")))
}

triplosens = function(all) {
    compg = all$gene[all$hit]
    ts = readxl::read_xlsx("../misc/triplosensitive_compare/1-s2.0-S0092867422007887-mmc7.xlsx")
    ts$is_comp = ifelse(ts$Gene %in% compg, "Compensated", "Other")
    wd = split(ts$pTriplo, ts$is_comp)
    wt = wilcox.test(wd[[1]], wd[[2]])

    ggplot(ts, aes(x=is_comp, y=pTriplo)) +
        geom_violin(aes(fill=is_comp), alpha=0.5) +
        scale_fill_manual(values=cm$cols[c("Compensated", "Other")], name="Genes") +
        ggbeeswarm::geom_quasirandom(data=ts[ts$is_comp=="Compensated",], alpha=0.5) +
        stat_summary(fun=mean, geom="crossbar", colour="red") +
        annotate("text", x=2.5, y=0.75, label=sprintf("p=%.2g (Wilcox)", wt$p.value)) +
        coord_flip() +
        theme_minimal() +
        theme(axis.text.y = element_blank()) +
        labs(title = "Triplosensitivity",
             x = "Compensation status",
             y = "Probability")
}

tcga_ccle_tissue = function() {
    cfg = yaml::read_yaml("../config.yaml")
    ccle = paste0(file.path("../model/fit_ccle-amp", cfg$ccle_tissues), ".rds") %>%
        lapply(readRDS) %>%
        setNames(cfg$ccle_tissues) %>% bind_rows(.id="tissue") %>%
        mutate(src = "CCLE")
    tcga = paste0(file.path("../model/fit_tcga_puradj-amp", cfg$tcga_tissues), ".rds") %>%
        lapply(readRDS) %>%
        setNames(cfg$tcga_tissues) %>% bind_rows(.id="tissue") %>%
        mutate(src = "TCGA")

    dset = bind_rows(ccle, tcga) %>%
        group_by(gene) %>% filter(all(c("CCLE", "TCGA") %in% src)) %>% ungroup() %>%
        mutate(tissue = ifelse(tissue == "pan", "Pan-Cancer", tissue),
               tissue = factor(tissue) %>% relevel("Pan-Cancer"),
               shrunk = estimate * (1-p.value),
               s = ifelse(shrunk < -0.3, 1, 0.7))

    res = bind_rows(dset, dset %>% mutate(src = "Common")) %>%
        mutate(src = factor(src, levels=c("Common", "CCLE", "TCGA"))) %>%
        group_by(src, gene) %>%
            filter(sum(!is.na(shrunk)) > 2) %>%
            summarize(n_comp = sum(shrunk < -0.3, na.rm=TRUE),
                      broom::tidy(lm(shrunk ~ 1)))
    sel = res %>% slice_max(n_comp, n=20, with_ties=FALSE) %>% select(gene, sel=src)

    dset2 = inner_join(dset, sel, relationship="many-to-many") %>%
        mutate(shrunk = pmax(-1, pmin(shrunk, 1)))
    p1 = ggplot(dset2, aes(x=gene, y=forcats::fct_rev(tissue), fill=shrunk)) +
        geom_tile(aes(width=s, height=s)) +
        scale_fill_distiller(palette="PuOr", name="Comp.\nscore") +
        facet_grid(src ~ sel, scales="free") +
        theme_minimal() +
        theme(strip.background = element_rect(color="black", linewidth=1),
              axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
        labs(x = "Gene",
             y = "Tissue")

    sets = gset$get_human("MSigDB_Hallmark_2020")
#    sres = split(res, res$src) %>%
#        lapply(gset$test_lm, sets=sets) %>%
#        bind_rows(.id="sel") %>% group_by(sel) %>% slice_min(adj.p, n=12, with_ties=FALSE) %>%
#        arrange(adj.p) %>% mutate(rank = factor(seq_len(n()), levels=seq_len(n()))) %>%
#        ungroup() %>% mutate(sel = factor(sel, levels=c("Common", "CCLE", "TCGA")))
    sres = res %>% filter(n_comp >= 3) %$% split(gene, src) %>%
        lapply(gset$test_fet, valid=unique(res$gene), sets=sets) %>%
        bind_rows(.id="sel") %>% group_by(sel) %>% slice_min(p.value, n=12, with_ties=FALSE) %>%
        arrange(p.value) %>% mutate(rank = factor(seq_len(n()), levels=seq_len(n()))) %>%
        ungroup() %>% mutate(sel = factor(sel, levels=c("Common", "CCLE", "TCGA")))

    p2 = ggplot(sres, aes(x=-log10(adj.p), y=forcats::fct_rev(rank))) +
        geom_col(fill="steelblue", alpha=0.2) +
        geom_text(aes(label=paste0(" ", label)), x=0, hjust=0) +
        facet_wrap(~ sel, scales="free") +
        theme_minimal() +
        theme(strip.background = element_rect(color="black", linewidth=1),
              axis.text.y = element_blank()) +
        labs(x = "-log10 FDR compensation in >= 3 Tissues",
             y = "Category")

    (p1 / p2) + plot_layout(heights=c(4,3))
}

sys$run({
    all = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")
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

    left = (tcga_vs_ccle() / go_cors()) + plot_layout(heights=c(1,3))
    right = (cna_comp(gistic, comp_all) / og_comp(comp) / rpe_comp(all) / triplosens(all)) +
        plot_layout(heights=c(1,1,1.5,0.8))

    asm = (((left | right) + plot_layout(widths=c(2,1))) / tcga_ccle_tissue()) +
        plot_layout(heights=c(6,5)) +
        plot_annotation(tag_levels='a') &
        theme(axis.text = element_text(size=10),
              legend.text = element_text(size=10),
              plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("FigS2-Compensation.pdf", 15, 20)
    print(asm)
    dev.off()
})
