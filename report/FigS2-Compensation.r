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
    ccle = readRDS("../model_compensation/fit_ccle-amp/pan.rds") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga1 = readRDS("../model_compensation/fit_tcga_naive-amp/pan.rds") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga2 = readRDS("../model_compensation/fit_tcga_pur-amp/pan.rds") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readRDS("../model_compensation/fit_tcga_puradj-amp/pan.rds") %>%
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
                   fill="#ffffffc0", label.size=NA, hjust=1, vjust=0.8) +
        ggnewscale::new_scale(c("fill")) +
        geom_point(data=hl, aes(fill=Gene), color="black", shape=21, size=2.5) +
        labs(y = "Compensation score TCGA",
             x = "Compensation score CCLE") +
        coord_cartesian(ylim=c(-1.2, 1.2)) +
        cm$theme_minimal() +
        theme(strip.background = element_rect(color=NA, fill="#ffffffc0"))
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
    use = grep("GO:00+(6260|6261|6302|8033|10467|6412|36297|381|43484|6397|48024|22618|6396)\\)$", both$label, value=TRUE)

    plt$denspt(both, aes(x=stat_tcga, y=stat_ccle, label=label),
               size=size_used, alpha=0.6, pal_alpha=0.5, tsize=4, n_tiles=80, h=40,
               draw_label=0, always_label=use) +
        scale_size_area(max_size=8, breaks=c(10,100,500,1000), name="Genes in set") +
        cm$theme_minimal() +
        guides(alpha="none") +
        labs(title = "Gene Ontology: Biological Process",
             x = "Compensation score TCGA (Wald stat.)",
             y = "Compensation score CCLE (Wald stat.)") +
        annotate("text", x=11, y=-6.5, color="blue", label=lab, parse=TRUE)
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
        labs(fill = "Frequent CNA", x = "Frequent CNA", y = "Δ ORF (Wald statistic)"),
        cm$theme_classic(),
        coord_cartesian(ylim=coordy, clip="off"),
        theme(axis.text.x = element_blank()),
        geom_hline(yintercept=median(y[both$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
    )

    p1 = ggplot(both, aes(x=type, y=estimate.x, fill=type)) +
        common(both$estimate.x, c(-0.3, 0.7), c(0.36, 0.55)) +
        labs(title = "CCLE", y="Compensation scores")
    p2 = ggplot(both, aes(x=type, y=estimate.y, fill=type)) +
        common(both$estimate.y, c(-0.5, 1.4), c(0.84, 1.2)) +
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
        labs(fill = "Driver status\n(freq. amplified)", x = "Driver status"),
        cm$theme_classic(),
        coord_cartesian(ylim=coordy, clip="off"),
        theme(axis.text.x = element_blank()),
        geom_hline(yintercept=median(y[both$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
    )

    p1 = ggplot(both, aes(x=type, y=estimate.x, fill=type)) +
        common(both$estimate.x, c(-0.3, 0.7), c(0.36, 0.55)) +
        labs(title = "CCLE", y="Compensation scores")
    p2 = ggplot(both, aes(x=type, y=estimate.y, fill=type)) +
        common(both$estimate.y, c(-0.5, 1.4), c(0.84, 1.2)) +
        labs(title = "TCGA", y="")

    (p1 | (p2 + plot_layout(tag_level="new"))) + plot_layout(guides="collect")
}

rpe_comp = function() {
    lookup = c(SS6="chr +7", SS51="+7 +22", SS111="+8 +9 +18")
    chrs = seq$gene_table() %>% select(label=external_gene_name, chr=chromosome_name) %>% distinct()
    comp = cm$get_comp_genes(pan=TRUE)
    dset = readRDS("../data/rnaseq_rpe1_broad/compute_fcs.rds") %>%
        tidyr::unnest(genes) %>%
        inner_join(chrs) %>%
        filter((term == "SS6" & chr == "7") |
               (term == "SS51" & chr %in% c("7", "22")) |
               (term == "SS111" & chr %in% c("8", "9", "18"))) %>%
        transmute(Sample = factor(lookup[term], levels=lookup),
                  Gene=label, chr=chr, LFC=log2FoldChange,
                  status = ifelse(Gene %in% comp, "Comp.", "Other"),
                  status = factor(status, levels=c("Other", "Comp."))) %>%
        group_by(Sample, chr) %>%
            mutate(LFC = scale(LFC, scale=FALSE)[,1]) %>%
        ungroup()

    ggplot(dset, aes(x=status, y=2^LFC, color=status)) +
        geom_boxplot(aes(fill=status), outlier.shape=NA, alpha=0.3) +
        ggbeeswarm::geom_quasirandom(dodge.width=0.8, aes(alpha=status)) +
        scale_y_log10() +
        facet_wrap(~ Sample) +
        coord_cartesian(ylim=c(0.5, 2)) +
        labs(x = "Isogenic RPE-1 clones",
             y = "Fold-change on\ngained chromosomes") +
        scale_color_manual(values=c(cm$cols[c("Other", "Comp.")]), name="Genes") +
        scale_fill_manual(values=c(cm$cols[c("Other", "Comp.")]), name="Genes") +
        scale_alpha_manual(values=c(Background=0.1, Compensated=0.6), guide="none") +
        cm$theme_minimal() +
        theme(axis.text.x = element_blank()) +
        ggsignif::geom_signif(color="black", y_position=-0.15,
            test=function(...) t.test(..., alternative="greater"),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Other", "Comp.")))
}

triplosens = function() {
    compg = cm$get_comp_genes(pan=TRUE)
    ts = readxl::read_xlsx("../misc/triplosensitive_compare/1-s2.0-S0092867422007887-mmc7.xlsx")
    ts$is_comp = ifelse(ts$Gene %in% compg, "Compensated", "Other")
    wd = split(ts$pTriplo, ts$is_comp)
    wt = wilcox.test(wd[[1]], wd[[2]])

    ggplot(ts, aes(x=is_comp, y=pTriplo)) +
        geom_violin(aes(fill=is_comp), color="#00000010", alpha=0.5) +
        scale_fill_manual(values=cm$cols[c("Compensated", "Other")], name="Genes") +
        ggbeeswarm::geom_quasirandom(data=ts[ts$is_comp=="Compensated",], alpha=0.5) +
        stat_summary(fun=mean, geom="crossbar", colour="red") +
        annotate("text", x=2.5, y=0.75, label=cm$fmt_p(wt$p.value), parse=TRUE) +
        coord_flip(clip="off") +
        cm$theme_minimal() +
        theme(axis.text.y = element_blank()) +
        labs(x = "Compensation\nstatus",
             y = "prob. Triplosensitivity")
}

venn_comp = function() {
    ov = readRDS("../misc/reviewer1/compensation.rds")$overlap
    names(ov)[names(ov) == "Goncalves"] = "Goncalves et al."
    names(ov)[names(ov) == "Schukken gene"] = "    Schukken et al. (gene)"
    names(ov)[names(ov) == "Schukken\nprotein"] = "Schukken et al.\n(protein)"
    ggvenn::ggvenn(ov, set_name_size=4, show_percentage=FALSE) +
        theme_void() + coord_cartesian(clip="off") +
        theme(axis.text.x = element_blank(), axis.text.y = element_blank())
}

tcga_mut = function(freqs) {
    ggplot(freqs, aes(x=class, y=freq)) +
        geom_boxplot(outlier.shape=NA) +
        ggbeeswarm::geom_quasirandom(aes(shape=class, fill=class),size=2, alpha=0.8) +
        scale_shape_manual(values=c(Other=NA, Compensated=21, ARGOS=21), name="Gene class") +
        scale_fill_manual(values=c(Other=NA, Compensated="#74ad9b", ARGOS="#de493d"), name="Gene class") +
        scale_y_log10() +
        labs(x = "Gene class",
             y = "TCGA mutation frequency") +
        ggsignif::geom_signif(color="black", y_position=c(-0.9,-0.6,-1), test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons = list(c("Other", "Compensated"),
                               c("Other", "ARGOS"),
                               c("Compensated", "ARGOS"))) +
        theme_classic() +
        coord_cartesian(clip="off")
}

rrm_pld = function() {
    res = readRDS("../misc/reviewer3/pld_domain.rds") %>%
        filter(label == "Compensated")
    ggplot(res, aes(x=estimate, y=p.value)) +
        geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), alpha=0.5, height=0.2) +
        geom_point(aes(shape=Comparison), fill="#74ad9b", size=3, alpha=0.8) +
        scale_x_log10() +
        scale_y_continuous(trans=ggforce::trans_reverser("log10"), breaks=c(0.05, 1e-10, 1e-20)) +
        scale_shape_manual(values=c(`RRM over all`=21, `PLD over RRM`=23), name="Compensated\ngenes") +
#        scale_fill_manual(guide=guide_legend(override.aes=list(shape=21)), name="Gene class",
#                          values=c(Toxic="#226b94", Compensated="#74ad9b", ARGOS="#de493d")) +
        geom_hline(yintercept=0.05, linetype="dashed") +
#        annotate("text", x=0.02, y=-log10(10), label=cm$fmt_p(0.05), vjust=-1, hjust=0, parse=TRUE) +
        labs(x = "Odds ratio (fold enrichment)",
             y = "P-value") +
        coord_cartesian(ylim=c(20,1e-27), clip="off") +
        cm$theme_classic() +
        theme(panel.grid.major = element_line(color="#dededea0"))
}

tcga_ccle_tissue = function() {
    dset = cm$get_comp_tissue() %>%
        mutate(s = ifelse(is_comp, 1, 0.7))

    res = bind_rows(dset, dset %>% mutate(src = "Common")) %>%
        mutate(src = factor(src, levels=c("Common", "CCLE", "TCGA"))) %>%
        group_by(src, gene) %>%
            filter(sum(!is.na(compensation)) > 2) %>%
            summarize(n_comp = sum(compensation < -0.3, na.rm=TRUE),
                      broom::tidy(lm(compensation ~ 1)))
    sel = res %>% slice_max(n_comp, n=20, with_ties=FALSE) %>% select(gene, sel=src)

    ov = readRDS("../misc/reviewer1/compensation.rds")$overlap
    names(ov)[names(ov) == "Goncalves"] = "Goncalves et al."
    names(ov)[names(ov) == "Schukken gene"] = "Schukken (gene)"
    names(ov)[names(ov) == "Schukken\nprotein"] = "Schukken (protein)"
    dset2 = inner_join(dset, sel, relationship="many-to-many") %>%
        mutate(src = paste(src, "data"),
               gene = factor(gene),
               compensation = pmax(-1, pmin(compensation, 1)))
    levels(dset2$sel) = paste(levels(dset2$sel), "top genes")
    ov2 = as_tibble(stack(ov)) %>%
        transmute(gene=factor(values, levels=levels(dset2$gene)), study=ind, in_study="Yes") %>%
        filter(grepl("Gonc|Schuk", study))
    dset3 = dset2 %>%
        select(gene, sel) %>% distinct() %>%
        inner_join(ov2, relationship="many-to-many") %>%
        mutate(tissue = droplevels(study), compensation = NA, src = "Study") %>%
        filter(!is.na(study))
    dset4 = bind_rows(dset2, dset3) %>%
        mutate(src = factor(src, levels=c("CCLE data", "TCGA data", "Study")))

    p1 = ggplot(dset4, aes(x=gene, y=forcats::fct_rev(tissue))) +
        geom_tile(aes(fill=compensation, width=s, height=s)) +
        scale_fill_distiller(palette="PuOr", name="Compensation\nscore", na.value="#ababab00") +
        ggnewscale::new_scale_fill() +
        geom_tile(aes(fill=in_study), color="white") +
        scale_fill_manual(values=c(Yes="#323232a0"), name="In Study",
                          na.value="transparent", na.translate=FALSE) +
        facet_grid(src ~ sel, scales="free", space="free") +
        cm$theme_minimal() +
        theme(strip.background = element_rect(color="black", linewidth=1),
              axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
        labs(x = "Gene",
             y = "Study / Tissue")

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
    levels(sres$sel) = paste(levels(sres$sel), "top genes")

    p2 = ggplot(sres, aes(x=-log10(adj.p), y=forcats::fct_rev(rank))) +
        geom_col(aes(fill=ifelse(estimate>1, "Enriched", "Depleted")), alpha=0.2) +
        scale_fill_manual(values=c(Enriched="steelblue", Depleted="coral"), name="Direction") +
        geom_text(aes(label=paste0(" ", label)), x=0, hjust=0) +
        facet_wrap(~ sel, scales="free") +
        cm$theme_minimal() +
        theme(strip.background = element_rect(color="black", linewidth=1),
              axis.text.y = element_blank()) +
        labs(x = "-log10 FDR compensation in ≥ 3 Tissues (Fisher's Exact Test)",
             y = "MSigDB Hallmark\ncategory")

    (p1/ p2) + plot_layout(heights=c(7,6))
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

    left = (tcga_vs_ccle() / go_cors()) + plot_layout(heights=c(1,3))
    right = (cna_comp(gistic, comp_all) / og_comp(comp) / rrm_pld() /
             wrap_elements(venn_comp() + theme(plot.margin=margin(0,-20,-10,-5,"mm")))) +
        plot_layout(heights=c(1,1,1.1,2.4))
    top = ((left | right) + plot_layout(widths=c(2,1)))

    mid = (tcga_mut(readRDS("../misc/reviewer3/mut_enrich.rds")) |
        rpe_comp() | triplosens()) + plot_layout(widths=c(5,6,4))

    asm = (top / mid / tcga_ccle_tissue()) +
        plot_layout(heights=c(6,1,5.5)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("FigS2-Compensation.pdf", 14.5, 19.5) # 14.5: extra space on labels left
    print(asm)
    dev.off()
})
