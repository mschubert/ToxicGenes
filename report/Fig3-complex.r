library(dplyr)
library(ggplot2)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
fig1 = import('./Fig1-Motivation')

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
    library("scales")
    reverselog_trans <- function(base = exp(1)) {
        trans <- function(x) -log(x, base)
        inv <- function(x) base^(-x)
        trans_new(paste0("reverselog-", format(base)), trans, inv,
                  log_breaks(base = base),
                  domain = c(1e-100, Inf))
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
    res$label[grepl("CPSF6|Cleavage", res$label)] = NA
    res$label = sub("components-", "", res$label)

    res2 = res %>%
        mutate(label = sub("(.*)", "`\\1`", label),
               label = ifelse(is.na(label) | !has_hit, label,
                              sprintf("%s^{%s}", label, hit_str)))

    ggplot(res2, aes(x=avg_orf, y=p.value)) +
        geom_rect(ymin=-Inf, ymax=2, xmin=-Inf, xmax=Inf, fill="#f3f3f3") +
        geom_rect(ymin=-Inf, ymax=Inf, xmin=-1, xmax=1, fill="#FAF4CD10") +
        geom_hline(yintercept=0.4, linetype="dashed", size=2, color="grey") +
        geom_vline(xintercept=0, linetype="dashed", size=2, color="grey") +
        geom_point(data=res %>% filter(!has_hit), aes(size=n, fill=has_hit), shape=21) +
        geom_point(data=res %>% filter(has_hit), aes(size=n, fill=has_hit), shape=21) +
        ggrepel::geom_label_repel(aes(label=label), max.overlaps=10, segment.alpha=0.3,
            label.size=NA, fill="#ffffffa0", min.segment.length=0, parse=TRUE,
            max.iter=1e5, max.time=10) +
        scale_fill_manual(values=c(`FALSE`="grey", `TRUE`="#FA524E")) +
        scale_size_binned_area(max_size=10) +
        scale_y_continuous(trans=reverselog_trans(10)) +
        xlim(c(max(res$avg_orf[res$p.value<0.2], na.rm=TRUE),
               min(res$avg_orf[res$p.value<0.2], na.rm=TRUE))) +
        theme_classic() +
        labs(x = "Mean ORF dropout compensated genes (Wald statistic)",
             y = "Overlap compensated genes (Fisher's Exact Test)")
}

sys$run({
    gistic_amp = fig1$get_gistic_scores() %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene=gene_name, frac)
    dset = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv") %>%
        inner_join(gistic_amp)

    boxes = wrap_elements(comp_hyp_box(dset))
    ov = overlap_venn(dset)
    cplx = complex_plot() #gistic_amp)

    asm = (((boxes / ov) + plot_layout(heights=c(1,1.5))) | cplx) +
        plot_layout(widths=c(1,1.8)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("Fig3-complex.pdf", 14, 8)
    print(asm)
    dev.off()
})
