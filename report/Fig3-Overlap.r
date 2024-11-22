library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
seq = import('seq')
gset = import('genesets')
cm = import('./common')

along_genome = function(dset, gistic, chrs=1:22) {
    lens = seq$chr_lengths("GRCh38", chrs=chrs) %>%
        stack() %>% transmute(x=1, xend=values, chr=ind) %>%
        rowwise() %>%
        mutate(scales = list(scale_x_continuous(limits=c(x, xend), expand=c(0,0))))

    cosmic = cm$get_cosmic_annot() %>% select(-tier, -hallmark) %>% filter(type != "OG+TSG")
    genes = gistic$genes %>% select(gene_name) %>% distinct() %>% mutate(type="Genes")
    comp = dset %>% filter(type == "Compensated") %>% select(type, gene_name=gene)
    hyp = dset %>% filter(type == "Hyperactivated") %>% select(type, gene_name=gene)
    orf = dset %>% filter(is_tox) %>% transmute(type="ORF dropout", gene_name=gene)

    labs = gistic$genes %>%
        filter(type == "amplification",
               gene_name %in% intersect(comp$gene_name, orf$gene_name)) %>%
        inner_join(gistic$smooth %>% select(type, chr, gam)) %>%
        rowwise() %>%
        mutate(frac = mgcv::predict.gam(gam, newdata=data.frame(tss=tss)))

    smooth = gistic$smooth %>% select(-gam) %>% tidyr::unnest(steps) %>%
        filter(type == "amplification", chr %in% chrs) %>%
        mutate(frac_amp = ifelse(frac > 0.15, frac, NA),
               type = stringr::str_to_title(type))
    sm_bg = smooth %>%
        group_by(chr) %>%
            mutate(seg_n = cumsum(is.na(frac_amp))) %>%
        na.omit() %>%
        group_by(chr, seg_n) %>%
            summarize(xmin=min(tss), xmax=max(tss)) %>%
        ungroup()

    dots = bind_rows(genes, cosmic, distinct(comp), distinct(hyp), orf) %>%
        mutate(type = factor(type, levels=unique(type))) %>%
        inner_join(gistic$genes %>% select(gene_name, chr, tss) %>% filter(!duplicated(gene_name))) %>%
        na.omit() %>%
        filter(chr %in% chrs)

    make_dens = function(tss, x, xend) {
        if (length(na.omit(tss)) == 0)
            return(tibble(x=c(x,xend), y=c(0,0)))
        dens = density(tss, from=unique(x), to=unique(xend), n=xend/1e6, bw=5e6,
            weights=rep(1, length(tss)), subdensity=TRUE)
        tibble(x=dens$x, y=dens$y)
    }
    res = tidyr::complete(dots, type, chr) %>%
        inner_join(lens) %>%
        group_by(type, chr) %>%
            arrange(tss) %>%
            summarize(dens = list(make_dens(tss, x, xend))) %>%
            tidyr::unnest(dens) %>%
        group_by(type) %>%
            mutate(y = y/max(y)) %>%
        ungroup()
    saveRDS(res, file="Fig3-Overlap.rds")
    labs2 = tidyr::expand_grid(labs %>% select(-type, -gam), distinct(res['type']))

    dens = ggplot(res, aes(x=x, y=y, fill=type)) +
        geom_rect(data=sm_bg, aes(xmin=xmin, xmax=xmax), ymin=-Inf, ymax=Inf, color=NA,
                  fill="firebrick", alpha=0.08, inherit.aes=FALSE) +
        geom_area(color="black", alpha=0.7) +
        geom_vline(data=labs2, aes(xintercept=tss), color="#656565", linetype="dashed", linewidth=0.5) +
        facet_grid(type ~ chr, scales="free", space="free") +
        theme_void() +
        theme(panel.background = element_rect(color=NA, fill="#efefef80"),
              panel.spacing.y = unit(0, "mm")) +
        scale_fill_manual(values=cm$cols[c("Genes", "Oncogene", "TSG",
            "Compensated", "Hyperactivated", "ORF dropout")], name="") +
        plot_layout(tag_level="new") +
        expand_limits(y=1.1) +
        coord_cartesian(clip="off", expand=FALSE) +
        guides(fill = guide_legend(byrow = TRUE)) +
        theme(strip.text.x = element_text(size=12),
              legend.title = element_text(size=12),
              legend.text = element_text(size=11),
              legend.spacing.y = unit(1.2, 'mm'))

    amp = ggplot(smooth, aes(x=tss)) +
        geom_segment(data=lens, aes(x=x, xend=xend), y=0, yend=0, alpha=1) +
        geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
        scale_fill_manual(values=cm$cols["Amplification"], name="CNA") +
        geom_line(aes(y=frac_amp, group=type, color="Frequently\namplified"),
                  lineend="round", size=1) +
        geom_point(data=labs, aes(x=tss, y=frac), color="black", fill="white", shape=21, size=2) +
        scale_color_manual(values=c("Frequently\namplified"="#960019"), name="") +
        facet_grid(. ~ chr, scales="free", space="free") +
        ggh4x::facetted_pos_scales(x=lens$scales) +
        theme_minimal() +
        guides(fill="none") +
        theme(legend.title = element_text(size=12),
              legend.text = element_text(size=11),
              axis.title = element_blank(),
              axis.text = element_blank(),
              strip.text.x = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        coord_cartesian(clip="off", expand=FALSE)

    labs3 = labs %>%
        mutate(tss2 = ifelse(gene_name == "SNRPA", tss+2e7, tss),
               gene_name = paste(gene_name, " "))
    glabs = ggplot(labs3, aes(x=tss, y=0)) +
        geom_segment(data=lens, aes(x=x, xend=xend), y=0.1, yend=0.1, alpha=0) +
        geom_segment(aes(xend=tss, yend=0.1)) +
        geom_text(aes(x=tss2, label=gene_name), size=3.5, angle=40, hjust=1, vjust=1) +
        facet_grid(. ~ chr, scales="free", space="free") +
        coord_cartesian(clip="off", ylim=c(-2,0.1), expand=FALSE) +
        theme_void() +
        theme(strip.text.x = element_blank(),
              panel.background = element_blank(),
              panel.spacing.y = unit(0, "mm")) +
    plot_layout(tag_level="new")

    (amp / dens / glabs) + plot_layout(heights=c(1,5,1.2)) &
        theme(plot.margin = unit(c(0,0,0,0), "mm"),
              strip.background = element_blank(),
              strip.text.y = element_blank(),
              panel.spacing.x = unit(1, "mm"))
}

comp_hyp_box = function(dset) {
    ds = dset %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = factor(type, levels=c("Background", "Compensated", "Hyperactivated")))

    ggplot(ds, aes(x=type, y=stat_orf, fill=type)) +
        geom_boxplot(outlier.shape=NA, alpha=0.7) +
        ggsignif::geom_signif(y_position=c(4.5, 6.5), color="black", test=t.test,
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            comparisons=list(c("Background", "Compensated"), c("Background", "Hyperactivated"))) +
        coord_cartesian(ylim=c(-7.5, 9)) +
        labs(fill = "Status", x = "Compensation status", y = "Δ ORF (Wald statistic)") +
        scale_fill_manual(values=cm$cols[c("Background", "Compensated", "Hyperactivated")]) +
        cm$theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(ds$stat_orf[ds$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
}

overlap_venn = function(dset, gistic_amp) {
    ov = list(CCLE = dset %>% filter(type_ccle == "Compensated") %>% pull(gene) %>% unique(),
              TCGA = dset %>% filter(type_tcga == "Compensated") %>% pull(gene) %>% unique(),
              ORF = dset %>% filter(is_tox) %>% pull(gene) %>% unique())
    all3 = Reduce(intersect, ov)
    all3 = ifelse(all3 %in% gistic_amp$gene, paste("▲", all3), all3)
    plt$venn(ov, alpha=0.4) +
        scale_fill_manual(values=cm$cols[c("TCGA", "CCLE", "ORF")]) +
        annotate("text", x=-13, y=13, label=paste(all3, collapse="\n"), size=4, hjust=1) +
        annotate("segment", x=-12.3, y=4.5, xend=-12.3, yend=21.5) +
        annotate("segment", x=-12, y=15, xend=4.5, yend=0.2) +
        coord_fixed(clip="off")
}

test_fet = function(set, corum, dset, comp, argos) {
    mat = matrix(nrow=2, c(
        length(intersect(corum[[set]], comp)),
        length(comp),
        length(corum[[set]]),
        length(dset$gene)
    ))
    mat[,2] = mat[,2] - mat[,1]
    broom::tidy(fisher.test(mat)) %>%
        mutate(n = length(corum[[set]]),
               avg_orf = mean(dset$stat_orf[dset$gene %in% corum[[set]]], na.rm=TRUE),
               hit_str = paste(sprintf("bold(`%s`)", intersect(argos, corum[[set]])), collapse=", "),
               has_hit = hit_str != "")
}

complex_plot = function(dset, hits) {
    .reverselog_trans = function(base=exp(1)) {
        scales::trans_new(paste0("log-", format(base)),
                          function(x) -log(x, base),
                          function(x) base^-x,
                          scales::log_breaks(base = base),
                          domain = c(1e-100, Inf))
    }
    .scientific_10 = function(x) {
        fmt = ifelse(x < 0.01, scales::scientific_format()(x), x)
        parse(text=gsub("1e", "10^", fmt))
    }

    corum = gset$get_human("CORUM_all") %>%
        gset$filter(min=3, valid=dset$gene)

    comp = dset %>% filter(is_comp) %>% pull(gene)
    argos = dset %>% filter(is_argos) %>% pull(gene)
    res = sapply(names(corum), test_fet, simplify=FALSE, corum=corum, dset=dset, comp=comp, argos=argos) %>%
        bind_rows(.id="set_name") %>%
        select(-method, -alternative) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) %>%
        mutate(label = ifelse(p.value < 0.002 | (has_hit & avg_orf < -1) |
                              (p.value < 0.1 & avg_orf < -4), set_name, NA))

    # names too long for nice alignment
    res$label[grepl("CPSF6|Cleavage|CDC5L|p130|SIN3B", res$label)] = NA
    res$label = sub("components-", "", res$label)
    res$p.value[res$label == "Large Drosha complex"] = 1e-9 # 1e-13

    res2 = res %>%
        mutate(label = sub("(.*)", "`\\1`", label),
               label = ifelse(is.na(label) | !has_hit, label,
                              sprintf("%s^{%s}", label, hit_str)))
    res2$label[res2$label == "`Large Drosha complex`"] = "`Large Drosha complex (`~italic(P)<10^{-10}~`)`"
    fdr = mean(c(res2$p.value[res$adj.p>0.2][1], rev(res2$p.value[res$adj.p<0.2])[1]))

    membs = corum[[grep("HEXIM1-DNA-PK", names(corum))]]
    gobp = gset$get_human("GO_Biological_Process_2021") %>% gset$filter(valid=membs)
    nhej = gobp[[grep("GO:0006303", names(gobp))]]
    ds2 = dset %>% filter(gene %in% membs) %>%
        select(gene, CCLE=comp_ccle, TCGA=comp_tcga, ORF=stat_orf) %>%
        mutate(gene = forcats::fct_reorder(gene, TCGA, .desc=TRUE)) %>%
        tidyr::gather("type", "value", -gene) %>%
        mutate(missing = ifelse(is.na(value), "No data", NA),
               type = factor(type, levels=c("TCGA", "CCLE", "ORF")),
               class = case_when(
                   gene %in% nhej ~ "NHEJ",
                   gene %in% c("SFPQ", "NONO", "PSPC1", "RBM14", "MATR3") ~ "Para-\nspeckle",
                   TRUE ~ "Other"),
               class = factor(class, levels=c("Para-\nspeckle", "Other", "NHEJ")))
    detail = ggplot(ds2, aes(x=value, y=gene, fill=class)) +
        geom_col() +
        geom_point(aes(shape=missing), x=0) +
        geom_vline(xintercept=0) +
        facet_wrap(~ type, scales="free_x") +
        theme_minimal() +
        theme(legend.key.size = unit(3, "mm"),
              legend.spacing.y = unit(-3, "mm"),
              axis.title.y = element_blank(),
              axis.title.x = element_text(size=10),
              legend.title = element_text(size=10),
              legend.text = element_text(size=10),
              strip.text.x = element_text(size=10),
              plot.background = element_rect(color="#e5e5e5", fill="#fdfdfd")) +
        scale_x_continuous(breaks=c(-0.5, -5)) +
        xlab("   Compensation (score) / ORF dropout (Wald)") +
        scale_fill_brewer(palette="Dark2", name="", direction=-1,
            guide=guide_legend(override.aes=list(shape=NA))) +
        scale_shape_manual(values=c("No data"=4), name="", na.translate=FALSE) +
        plot_layout(tag_level="new")

    assocs = ggplot(res2, aes(x=avg_orf, y=p.value)) +
        geom_rect(ymin=-Inf, ymax=1, xmin=-Inf, xmax=Inf, fill="#f3f3f3") +
        geom_rect(ymin=-Inf, ymax=Inf, xmin=-1, xmax=1, fill="#FAF4CD10") +
        geom_vline(xintercept=0, linetype="dashed", size=2, color="grey") +
        geom_hline(yintercept=fdr, linetype="dashed", color="black") +
        annotate("text", y=fdr, x=-5.25, vjust=-1, hjust=0, label="20% FDR", size=3.5) +
        geom_point(data=res %>% filter(!has_hit), aes(size=n, fill=has_hit), shape=21) +
        geom_point(data=res %>% filter(has_hit), aes(size=n, fill=has_hit), shape=21) +
        ggrepel::geom_label_repel(aes(label=label), max.overlaps=12, segment.alpha=0.3,
            label.size=NA, fill="#ffffffd0", min.segment.length=0, parse=TRUE,
            max.iter=1e5, max.time=10, seed=1) +
        scale_fill_manual(values=c(`FALSE`="grey", `TRUE`=cm$cols[["Comp+ORF"]])) +
        guides(fill = guide_legend(override.aes=list(size=3))) +
        scale_size_binned_area(max_size=10) +
        scale_y_continuous(trans=.reverselog_trans(10), labels=.scientific_10) +
        xlim(c(max(res$avg_orf[res$p.value<0.2], na.rm=TRUE),
               min(res$avg_orf[res$p.value<0.2], na.rm=TRUE))) +
        cm$theme_classic() +
        labs(x = "Mean ORF dropout compensated genes (Wald statistic)",
             y = "Overlap compensated genes (p-value Fisher's Exact Test)",
             size = "Protein\ncomplex\nmembers",
             fill = "Contains\nARGOS\ngene")

    assocs +
        annotate("curve", x=-1.8, y=1.5e-5, xend=-2.15, yend=4e-6, color="black",
                 curvature=-0.4, lineend="round", linejoin="round",
                 arrow=arrow(type="closed", length=unit(2.5,"mm"))) +
        inset_element(detail, 0.54, 0.57, 1, 0.84)
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")
    gistic_amp = gistic$genes %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene=gene_name, frac)
    dset = cm$get_pancan_summary() %>% left_join(gistic_amp)

    top = along_genome(dset, gistic)
    boxes = wrap_elements(comp_hyp_box(dset))
    ov = wrap_elements(overlap_venn(dset, gistic_amp) + theme(plot.margin = unit(c(0,-15,-10,0), "mm")))
    cplx = complex_plot(dset)

    asm = (wrap_plots(top) /
        ((((boxes / ov) + plot_layout(heights=c(1,1.5))) | cplx) +
        plot_layout(widths=c(1,1.8))) + plot_layout(heights=c(1,2.3))) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=24, face="bold"))

    cairo_pdf("Fig3-Overlap.pdf", 14, 11)
    print(asm)
    dev.off()
})
