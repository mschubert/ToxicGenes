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
        labs(title = "CCLE", y="Compensation\nscores")
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
        labs(title = "CCLE", y="Compensation\nscores")
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

    ggplot(dset, aes(x=status, y=2^LFC, fill=status)) +
        geom_boxplot(color="#40404055", outlier.shape=NA) +
        ggbeeswarm::geom_quasirandom(data=dset[dset$status == "Comp.",],
            shape=21, color="black", alpha=0.9, dodge.width=0.8) +
        scale_y_log10() +
        facet_wrap(~ Sample) +
        coord_cartesian(ylim=c(0.5, 2)) +
        labs(x = "Isogenic RPE-1 clones",
             y = "Fold-change on\ngained chromosomes") +
        scale_fill_manual(values=cm$col_study, guide="none") +
        cm$theme_minimal() +
        theme(axis.text.x = element_blank()) +
        ggsignif::geom_signif(color="black", y_position=-0.15,
            test=function(...) t.test(..., alternative="greater"),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Other", "Comp.")))
}

rpe_study = function() {
    chrs = seq$gene_table() %>% select(label=external_gene_name, chr=chromosome_name) %>% distinct()
    lookup = c(all="All genes", ours="ours", Goncalves="Goncalves",
               `Schukken\nprotein`="Schukken\n(protein)", `Schukken gene`="Schukken\n(gene)")
    ov = stack(readRDS("../misc/reviewer1/compensation.rds")$overlap) %>%
        dplyr::rename(Gene=values, Study=ind) %>%
        mutate(Study = lookup[as.character(Study)])
    dset1 = readRDS("../data/rnaseq_rpe1_broad/compute_fcs.rds") %>%
        tidyr::unnest(genes) %>%
        inner_join(chrs) %>%
        filter((term == "SS6" & chr == "7") |
               (term == "SS51" & chr %in% c("7", "22")) |
               (term == "SS111" & chr %in% c("8", "9", "18"))) %>%
        transmute(Sample=term, Gene=label, chr=chr, LFC=log2FoldChange) %>%
        group_by(Sample, chr) %>%
            mutate(LFC = scale(LFC, scale=FALSE)[,1]) %>%
        ungroup()
    dset2 = dset1 %>% inner_join(ov, relationship="many-to-many") %>%
        bind_rows(dset1 %>% mutate(Study = "All genes")) %>%
        mutate(group = ifelse(grepl("Gonc|Schukk", Study), "prev", Study),
               Study = factor(Study, levels=lookup))

    ggplot(dset2, aes(x=group, fill=Study, y=2^LFC)) +
        geom_boxplot(outlier.shape=NA) +
        scale_y_log10() +
        scale_fill_manual(values=cm$col_study) +
        coord_cartesian(ylim=c(0.3, 4), clip="off") +
        labs(x = "Study",
             y = "Fold-change on\ngained chromosomes") +
        ggsignif::geom_signif(color="black", y_position=c(0.48, 0.12, 0.28),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0, test=t.test,
            comparisons = list(c("All genes", "prev"),
                               c("prev", "ours"),
                               c("All genes", "ours"))) +
        cm$theme_classic() +
        theme(axis.text.x = element_blank())
}

triplosens = function() {
    ts = readxl::read_xlsx("../misc/triplosensitive_compare/1-s2.0-S0092867422007887-mmc7.xlsx") %>%
        select(gene=Gene, pTriplo)
    lookup = c(all="All genes", ours="ours", Goncalves="Goncalves",
               `Schukken\nprotein`="Schukken\n(protein)", `Schukken gene`="Schukken\n(gene)")
    ov = stack(readRDS("../misc/reviewer1/compensation.rds")$overlap) %>%
        transmute(gene=values, Study=lookup[as.character(ind)])
    dset = inner_join(ov, ts) %>%
        bind_rows(ts %>% mutate(Study="All genes")) %>%
        mutate(group = ifelse(grepl("Gonc|Schukk", Study), "prev", Study),
               Study = factor(Study, levels=lookup))
    ggplot(dset, aes(x=group, fill=Study, y=pTriplo)) +
        geom_boxplot() +
        labs(x = "Study",
             y = "Probability of\nTriplosensitivity") +
        scale_y_continuous(breaks=c(0, 0.5, 1)) +
        scale_fill_manual(values=cm$col_study) +
        coord_cartesian(ylim=c(0,1.15), clip="off") +
        cm$theme_classic() +
        theme(axis.text.x = element_blank()) +
        ggsignif::geom_signif(color="black", y_position=c(1.4, 1.03, 1.19),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0, test=t.test,
            comparisons = list(c("All genes", "prev"),
                               c("prev", "ours"),
                               c("All genes", "ours")))
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

tcga_mut = function() {
    lookup = c(all="All genes", ours="ours", Goncalves="Goncalves",
               `Schukken\nprotein`="Schukken\n(protein)", `Schukken gene`="Schukken\n(gene)")
    ov = stack(readRDS("../misc/reviewer1/compensation.rds")$overlap) %>%
        dplyr::rename(gene=values, coll=ind)
    freqs = readRDS("../misc/reviewer3/mut_enrich.rds") %>% select(-class)
    freqs2 = left_join(freqs, ov %>% bind_rows(data.frame(gene=freqs$gene, coll="all"))) %>%
        mutate(group = ifelse(grepl("Gonc|Schukk", coll), "prev", coll),
               coll = factor(lookup[coll], levels=lookup))

    ggplot(freqs2, aes(x=group, fill=coll, y=freq)) +
        geom_boxplot(outlier.shape=NA) +
        scale_y_log10() +
        scale_fill_manual(values=cm$col_study) +
        ggsignif::geom_signif(color="black", y_position=c(-0.6,-1.25,-0.97), test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons = list(c("all", "prev"),
                               c("prev", "ours"),
                               c("all", "ours"))) +
        labs(x = "Study",
             y = "TCGA mutation\nfrequency",
             fill = "Study") +
        cm$theme_classic() +
        theme(axis.text.x = element_blank()) +
        coord_cartesian(ylim=c(1e-3, 0.1), clip="off")
}

rrm_pld = function() {
    gs = readRDS("../misc/reviewer3/pld_domain.rds")[c("rrm", "prion")]
    comp = list(comp=cm$get_comp_genes(pan=TRUE))
    gclass = seq$gene_table() %>% select(external_gene_name, gene_biotype) %>% distinct() %>%
        filter(gene_biotype %in% c("lncRNA", "processed_pseudogene", "protein_coding"),
               !is.na(external_gene_name)) %>%
        unstack()

    res = bind_rows(.id="comp",
        `protein coding` = gset$test_fet(valid=comp_all$gene_name, hits=gclass$protein_coding, sets=comp),
        lncRNA = gset$test_fet(valid=comp_all$gene_name, hits=gclass$lncRNA, sets=comp),
        pseudogene =  gset$test_fet(valid=comp_all$gene_name, hits=gclass$processed_pseudogene, sets=comp),
        `RRM vs.\nprotein coding` =  gset$test_fet(valid=gclass$protein_coding, hits=gs$rrm, sets=comp),
        `PLD vs. RRM` =  gset$test_fet(valid=gs$rrm, hits=gs$prion, sets=comp)
    ) %>% mutate(comp = factor(comp, levels=comp))

    ggplot(res, aes(x=estimate, y=p.value)) +
        geom_vline(xintercept=1, color="darkgrey", linetype="dotted") +
        geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), alpha=0.5, height=0.2) +
        geom_point(aes(color=comp), fill="#74ad9b", size=3, alpha=0.8) +
        scale_x_log10(breaks=c(`0.001`=0.001, `0.01`=0.01, `0.1`=0.1, `1`=1, `10`=10, `30`=30)) +
        scale_y_continuous(trans=ggforce::trans_reverser("log10"), breaks=c(0.05, 1e-10, 1e-20, 1e-30)) +
        geom_hline(yintercept=0.05, linetype="dashed") +
        labs(x = "Odds ratio (fold enrichment)",
             y = "P-value",
             color = "Comparison") +
        coord_cartesian(ylim=c(20,NA), clip="off") +
        cm$theme_classic() +
        theme(panel.grid.major = element_line(color="#dededea0"))
}

tcga_ccle_tissue = function() {
    dset = cm$get_comp_tissue() %>%
        mutate(s = ifelse(!is.na(type) & type == "Compensated", 1, 0.7))

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
    sres = res %>% filter(n_comp >= 3) %$% split(gene, src) %>%
#        c(ov[c("Goncalves et al.", "Schukken (protein)", "Schukken (gene)")]) %>%
        lapply(gset$test_fet, valid=unique(res$gene), sets=sets)
    sres2 = sres[c("Common", "CCLE", "TCGA")] %>% bind_rows(.id="sel") %>%
        group_by(sel) %>% slice_min(p.value, n=12, with_ties=FALSE) %>%
        arrange(p.value) %>% mutate(rank = factor(seq_len(n()), levels=seq_len(n()))) %>%
        ungroup() %>% mutate(sel = factor(sel, levels=c("Common", "CCLE", "TCGA")))
    levels(sres2$sel) = paste(levels(sres2$sel), "top genes")
#    add = sres[c("Goncalves et al.", "Schukken (protein)", "Schukken (gene)")] %>%
#        bind_rows(.id = "Study\ncomparison") %>%
#        mutate(sel = factor("Common top genes"),
#               rank = factor(match(label, sres2$label[sres2$sel == "Common top genes"])),
#               adj.p = pmax(adj.p, 1e-10)) %>%
#        filter(label %in% sres2$label[sres2$sel == "Common top genes"])

#    p2 = ggplot(sres2, aes(x=abs(log2(estimate)), y=forcats::fct_rev(rank))) +
    p2 = ggplot(sres2, aes(x=-log10(p.value), y=forcats::fct_rev(rank))) +
        geom_col(aes(fill=ifelse(estimate>1, "Enriched", "Depleted")), alpha=0.2) +
        scale_fill_manual(values=c(Enriched="steelblue", Depleted="coral"), name="Direction") +
#        geom_point(data=add, aes(color=`Study\ncomparison`), fill=NA, size=3, alpha=0.7) +
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

    ccle = readxl::read_xlsx("SuppData1_CCLE-comp.xlsx", sheet="Pan-Cancer") %>%
        mutate(estimate = pmax(-2, pmin(compensation, 2.5)))
    tcga3 = readxl::read_xlsx("SuppData2_TCGA-comp.xlsx", sheet="Pan-Cancer") %>%
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

    mid = (tcga_mut() | rpe_comp() | rpe_study() | triplosens()) +
        plot_layout(widths=c(2,3,2,2), guides="collect")

    asm = (top / mid / tcga_ccle_tissue()) +
        plot_layout(heights=c(6,1,5.5)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("FigS2-Compensation.pdf", 14.5, 19.5) # 14.5: extra space on labels left
    print(asm)
    dev.off()
})
