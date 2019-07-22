library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'rds', 'merge_genes.rds'),
    opt('y', 'yaml', 'yaml', 'top_genes.yaml'),
    opt('p', 'plotfile', 'pdf', 'top_genes.pdf'))

#' Get the percentile of x in y
plot_stats = function(gene) {
    cur = dset %>%
        group_by(dset, fit, adj) %>%
        mutate(pctile = 100 * (1-rank(statistic)/n())) %>%
        filter(name == gene)

    ggplot(dset, aes(x=1, y = statistic, color=adj)) +
        geom_hline(yintercept=0, linetype="dashed", color="grey") +
        geom_violin(position="identity", alpha=0) +
        geom_point(data=cur, size=5, alpha=0.7) +
        ggrepel::geom_text_repel(data=cur, size=2, box.padding=unit(7, "pt"),
            aes(label=sprintf("%.2f th\np %.1g", pctile, adj.p)),
            parse=FALSE, color="black") +
        facet_wrap(~ dset + fit, scale="free_x", nrow=1) +
        coord_cartesian(ylim=c(0,-40)) +
        guides(color = FALSE) +
        labs(title = gene) +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
}

dset = readRDS(args$dset) %>%
    mutate(dset = relevel(factor(dset), "orf"))
top = yaml::read_yaml(args$yaml)$genes #TODO: use right set if not only genes

overview = lapply(top, plot_stats)

pdf(args$plotfile, 14, 12)
cowplot::plot_grid(plotlist=overview)
dev.off()
