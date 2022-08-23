library(dplyr)
library(ggplot2)
sys = import('sys')
gset = import('genesets')
fig1 = import('./Fig1-Motivation')

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
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
               has_hit = length(intersect(hits, corum[[set]])) != 0)
}

complex_plot = function(gistic_amp,
        label_sets=c("SMN complex", "PCNA-p21 complex", "OCT1-SOX2 DNA-protein complex")) {
    dset = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv") %>%
        inner_join(gistic_amp)
    corum = gset$get_human("CORUM_all") %>%
        gset$filter(min=1, valid=dset$gene)

    res = sapply(names(corum), test_fet, simplify=FALSE, corum=corum, dset=dset) %>%
        bind_rows(.id="set_name") %>%
        select(-method, -alternative) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               label = ifelse(p.value < 5e-4 | set_name %in% label_sets, set_name, NA)) %>%
        arrange(adj.p, p.value)

    ggplot(res, aes(x=avg_orf, y=p.value)) +
        geom_rect(ymin=-Inf, ymax=2, xmin=-Inf, xmax=Inf, fill="#f3f3f3") +
        geom_rect(ymin=-Inf, ymax=Inf, xmin=-1, xmax=1, fill="#FAF4CD10") +
        geom_hline(yintercept=0.05, linetype="dashed", size=2, color="grey") +
        geom_vline(xintercept=0, linetype="dashed", size=2, color="grey") +
        geom_point(aes(size=n, fill=has_hit), shape=21) +
        ggrepel::geom_text_repel(aes(label=label)) +
        scale_fill_manual(values=c(`FALSE`="grey", `TRUE`="#FA524E")) +
        scale_size_binned_area(max_size=10) +
        scale_y_continuous(trans=reverselog_trans(10)) +
        scale_x_reverse() +
        theme_classic() +
        labs(x = "Mean ORF dropout compensated genes (Wald statistic)",
             y = "Overlap compensated genes (Fisher's Exact Test)")
}

sys$run({
    gistic_amp = fig1$get_gistic_scores() %>%
        filter(type == "amplification", frac > 0.15) %>%
        select(gene=gene_name, frac)

    p = complex_plot(gistic_amp)

    pdf("Fig3-complex.pdf", 8, 6)
    print(p)
    dev.off()
})
