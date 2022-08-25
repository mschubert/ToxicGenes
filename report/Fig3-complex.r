library(dplyr)
library(ggplot2)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
fig1 = import('./Fig1-Motivation')

along_genome = function(gistic) {
    cosmic = fig1$get_cosmic_annot() %>% select(-tier)
    genes = gistic$genes %>% select(gene_name) %>% distinct() %>% mutate(type="Genes")

    dset = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")
    comp = tibble(type="Compensated", gene_name=dset$gene[dset$est_ccle < -0.3 & dset$est_tcga < -0.3])
    hyp = tibble(type="Hyperactivated", gene_name=dset$gene[dset$est_ccle > 0.3 & dset$est_tcga > 0.3])
    orf = tibble(type="ORF dropout", gene_name=dset$gene[dset$stat_orf < -5])

    smooth = gistic$smooth %>% select(-gam) %>% tidyr::unnest(steps) %>%
        filter(type == "amplification") %>%
        mutate(frac_amp = ifelse(frac > 0.15, frac, NA),
               type = stringr::str_to_title(type))

    dots = bind_rows(genes, cosmic, comp, hyp, orf) %>%
        mutate(type = factor(type, levels=unique(type))) %>%
        inner_join(gistic$genes %>% select(gene_name, chr, tss) %>% distinct()) %>%
        na.omit()
    dens = ggplot(dots, aes(x=tss, y=type, fill=type)) +
        ggridges::geom_density_ridges(scale=0.9, bandwidth=5e6, jittered_points=TRUE,
            position = ggridges::position_points_jitter(width=0.05, height=0),
            point_shape='|', point_size=2, point_alpha=1, alpha=0.7) +
        facet_grid(. ~ chr, scales="free", space="free") +
        theme_void() +
        theme(panel.spacing.x = unit(1, "mm")) +
        coord_cartesian(clip="off", expand=FALSE) +
        scale_y_discrete(limits=rev) +
        scale_fill_discrete(name="") +
        plot_layout(tag_level="new")

    amp = ggplot(smooth, aes(x=tss)) +
        geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
        scale_fill_manual(values=c(Amplification="firebrick"), name="CNA") +
        geom_line(aes(y=frac_amp, group=type, color="Frequently\namplified"),
                  lineend="round", size=1) +
        scale_color_manual(values=c("Frequently\namplified"="#960019"), name="") +
        facet_grid(. ~ chr, scales="free", space="free") +
        theme_minimal() +
        guides(fill="none") +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.spacing.x = unit(1, "mm")) +
        coord_cartesian(clip="off", expand=FALSE)

    (amp / dens) + plot_layout(heights=c(1,4))
}

comp_hyp_box = function(dset) {
    ds = dset %>%
        mutate(type = case_when(
            est_ccle < -0.3 & est_tcga < -0.3 ~ "Compensated",
            est_ccle > 0.3 & est_tcga > 0.3 ~ "Hyperactivated",
            TRUE ~ "Background"
        )) %>%
        mutate(type = factor(type, levels=c("Background", "Compensated", "Hyperactivated")))

    ggplot(ds, aes(x=type, y=stat_orf, fill=type)) +
        geom_boxplot(outlier.shape=NA) +
        ggsignif::geom_signif(y_position=c(4.5, 6.5), color="black", test=wilcox.test,
            comparisons=list(c("Background", "Compensated"), c("Background", "Hyperactivated"))) +
        coord_cartesian(ylim=c(-7, 9)) +
        labs(fill = "Status", x = "Compensation set", y = "Δ ORF (Wald statistic)") +
        theme_classic() +
        theme(axis.text.x = element_blank()) +
        geom_hline(yintercept=median(ds$stat_orf[ds$type=="Background"], na.rm=TRUE),
                   linetype="dashed", color="black")
}

overlap_venn = function(dset) {
    ov = list(CCLE = unique(dset$gene[dset$est_ccle < -0.3]),
              TCGA = unique(dset$gene[dset$est_tcga < -0.3]),
              ORF = unique(dset$gene[dset$stat_orf < -5 & !is.na(dset$stat_orf)]))
    all3 = Reduce(intersect, ov)
    plt$venn(ov, alpha=0.4) +
        scale_fill_manual(values=c(TCGA="#74ad9b", CCLE="#226b94", ORF="#f7974e")) +
        annotate("text", x=-8, y=8, label=paste(all3, collapse="\n"), size=4, hjust=1) +
        annotate("segment", x=-7.5, y=4, xend=-7.5, yend=11.8) +
        annotate("segment", x=-7, y=8, xend=2.8, yend=-0.1)
}

test_fet = function(set, corum, dset, hits=c("RBM14", "POU2F1", "CDKN1A", "SNRPA", "ZBTB14")) {
    mat = matrix(nrow=2, c(
        length(intersect(corum[[set]], dset$gene[dset$hit])),
        length(dset$gene[dset$hit]),
        length(corum[[set]]),
        length(dset$gene)
    ))
    mat[,2] = mat[,2] - mat[,1]
    broom::tidy(fisher.test(mat)) %>%
        mutate(n = length(corum[[set]]),
               avg_orf = mean(dset$stat_orf[dset$gene %in% corum[[set]]], na.rm=TRUE),
               hit_str = paste(sprintf("bold(`%s`)", intersect(hits, corum[[set]])), collapse=", "),
               has_hit = hit_str != "")
}

complex_plot = function() {
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

    dset = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")
    corum = gset$get_human("CORUM_all") %>%
        gset$filter(min=3, valid=dset$gene)

    res = sapply(names(corum), test_fet, simplify=FALSE, corum=corum, dset=dset) %>%
        bind_rows(.id="set_name") %>%
        select(-method, -alternative) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) %>%
        mutate(label = ifelse(p.value < 0.002 | (has_hit & avg_orf < -1) |
                              (p.value < 0.1 & avg_orf < -4), set_name, NA))

    # names too long for nice alignment
    res$label[grepl("CPSF6|Cleavage|CDC5L", res$label)] = NA
    res$label = sub("components-", "", res$label)
    res$p.value[res$label == "Large Drosha complex"] = 1e-9 # 1e-13

    res2 = res %>%
        mutate(label = sub("(.*)", "`\\1`", label),
               label = ifelse(is.na(label) | !has_hit, label,
                              sprintf("%s^{%s}", label, hit_str)))
    res2$label[res2$label == "`Large Drosha complex`"] = "`Large Drosha complex (`~p<10^{-10}~`)`"
    fdr = mean(c(res2$p.value[res$adj.p>0.2][1], rev(res2$p.value[res$adj.p<0.2])[1]))

    ggplot(res2, aes(x=avg_orf, y=p.value)) +
        geom_rect(ymin=-Inf, ymax=1, xmin=-Inf, xmax=Inf, fill="#f3f3f3") +
        geom_rect(ymin=-Inf, ymax=Inf, xmin=-1, xmax=1, fill="#FAF4CD10") +
        geom_vline(xintercept=0, linetype="dashed", size=2, color="grey") +
        geom_hline(yintercept=fdr, linetype="dashed", color="black") +
        annotate("text", y=fdr, x=-5.45, vjust=-1, hjust=0, label="20% FDR", size=3) +
        geom_point(data=res %>% filter(!has_hit), aes(size=n, fill=has_hit), shape=21) +
        geom_point(data=res %>% filter(has_hit), aes(size=n, fill=has_hit), shape=21) +
        ggrepel::geom_label_repel(aes(label=label), max.overlaps=10, segment.alpha=0.3,
            label.size=NA, fill="#ffffffa0", min.segment.length=0, parse=TRUE,
            max.iter=1e5, max.time=10) +
        scale_fill_manual(values=c(`FALSE`="grey", `TRUE`="#FA524E")) +
        scale_size_binned_area(max_size=10) +
        scale_y_continuous(trans=.reverselog_trans(10), labels=.scientific_10) +
        xlim(c(max(res$avg_orf[res$p.value<0.2], na.rm=TRUE),
               min(res$avg_orf[res$p.value<0.2], na.rm=TRUE))) +
        theme_classic() +
        labs(x = "Mean ORF dropout compensated genes (Wald statistic)",
             y = "Overlap compensated genes (p-value Fisher's Exact Test)",
             size = "Protein\ncomplex\nmembers",
             fill = "Contains\ngene of\ninterest")
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")
    gistic_amp = gistic$genes %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene=gene_name, frac)
    dset = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv") %>%
        inner_join(gistic_amp)

    top = along_genome(gistic)
    boxes = wrap_elements(comp_hyp_box(dset))
    ov = overlap_venn(dset)
    cplx = complex_plot()

    asm = (wrap_plots(top) /
        ((((boxes / ov) + plot_layout(heights=c(1,1.5))) | cplx) +
        plot_layout(widths=c(1,1.8))) + plot_layout(heights=c(1,3.4))) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("Fig3-complex.pdf", 14, 11)
    print(asm)
    dev.off()
})
