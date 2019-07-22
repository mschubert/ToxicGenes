library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'rds', 'merge_genes.rds'),
    opt('y', 'yaml', 'yaml', 'top_genes.yaml'),
    opt('p', 'plotfile', 'pdf', 'top_genes.pdf'))

do_plot = function(gene) {
    dset %>%
        filter(name == gene) %>%
        mutate(dset = relevel(factor(dset), "orf")) %>%
        ggplot(aes(x=name, y = statistic, color=adj)) +
            geom_point(size=5) +
            facet_wrap(~ dset + fit, scale="free_x", nrow=1) +
            expand_limits(y=0) +
            geom_hline(yintercept=0, linetype="dashed", color="grey") +
            guides(color = FALSE) +
            labs(title = gene) +
            theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank())
}

dset = readRDS(args$dset)
top = yaml::read_yaml(args$yaml)$genes #TODO: use right set if not only genes

plots = lapply(top, do_plot)

pdf(args$plotfile, 12, 10)
cowplot::plot_grid(plotlist=plots)
dev.off()
