library(dplyr)
library(ggplot2)
sys = import('sys')
gset = import('genesets')

library("scales")
reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

test_fet = function(set) {
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
               has_RBM14 = "RBM14" %in% corum[[set]])
}

complex_plot = function() {
    dset = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")
    corum = gset$get_human("CORUM_all") %>%
        gset$filter(min=5)

    res = sapply(names(corum), test_fet, simplify=FALSE) %>%
        bind_rows(.id="set_name") %>%
        select(-method, -alternative) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)

    ggplot(res, aes(x=avg_orf, y=p.value)) +
        geom_rect(ymin=-Inf, ymax=2, xmin=-Inf, xmax=Inf, fill="#f3f3f3") +
        geom_rect(ymin=-Inf, ymax=Inf, xmin=-1, xmax=1, fill="#FAF4CD10") +
        geom_hline(yintercept=0.05, linetype="dashed", size=2, color="grey") +
        geom_vline(xintercept=0, linetype="dashed", size=2, color="grey") +
        geom_point(aes(size=n, fill=has_RBM14), shape=21) +
        ggrepel::geom_text_repel(aes(label=set_name), max.overlaps=5) +
        scale_fill_manual(values=c(`FALSE`="grey", `TRUE`="#FA524E")) +
        scale_size_binned_area(max_size=10) +
        scale_y_continuous(trans=reverselog_trans(10)) +
        scale_x_reverse() +
        theme_classic() +
        labs(x = "Mean ORF dropout compensated genes (Wald statistic)",
             y = "Overlap compensated genes (Fisher's Exact Test)")
}

sys$run({
    p = complex_plot()

    pdf("Fig4-complex.pdf", 8, 6)
    print(p)
    dev.off()
})
