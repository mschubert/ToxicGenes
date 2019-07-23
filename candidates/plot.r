library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'rds', 'merge/pan/genes.rds'),
    opt('y', 'yaml', 'yaml', 'pan/top-genes.yaml'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('p', 'plotfile', 'pdf', 'pan/top-genes.pdf'))

#' Get the percentile of x in y
plot_stats = function(gene) {
    cur = filter(dset, name == gene)
    ggplot(dset, aes(x=1, y = statistic, color=adj)) +
        geom_hline(yintercept=0, linetype="dashed", color="grey") +
        geom_violin(position="identity", alpha=0) +
        geom_point(data=cur, size=5, alpha=0.7) +
        ggrepel::geom_text_repel(data=cur, size=2, box.padding=unit(7, "pt"),
            aes(label=sprintf("%.2f th\nFDR %.1g", pctile, adj.p)),
            parse=FALSE, color="black", direction="y", segment.alpha=0.3) +
        facet_wrap(~ dset + fit, scale="free_x", nrow=1) +
        coord_cartesian(ylim=c(0,min(cur$statistic, na.rm=TRUE))) +
        guides(color = FALSE) +
        labs(title = gene) +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
}

top = yaml::read_yaml(args$yaml)$genes #TODO: use right set if not only genes

#TODO: get the min of all top hits so we set limits?

dset = readRDS(args$dset)
overview = lapply(top, plot_stats)

###
### ORF data
###
orfdata = readRDS("../orf/overview.rds")
if (args$tissue != "pan")
    orfdata = filter(orfdata, tissue == args$tissue)
orfdata = orfdata %>%
    mutate(`GENE SYMBOL` = factor(`GENE SYMBOL`, levels=top)) %>%
    group_by(`GENE SYMBOL`) %>%
        mutate(construct_i = as.integer(factor(`Construct IDs`))) %>%
    ungroup()
porf = ggplot(orfdata, aes(x=DMSO, y=z_LFC)) +
    geom_hline(yintercept=0, color="red") +
    geom_hline(yintercept=c(-1,1), color="red", linetype="dotted") +
    geom_point(aes(color=tissue, shape=factor(construct_i)), size=3, alpha=0.6) +
    facet_wrap(~ `GENE SYMBOL`) +
    ggtitle("ORF drop-out (loess normalized, red line: mean +/- SD)")

###
### CCLE data
###
ccledata = readRDS("../data/ccle/dset.rds")
names(dimnames(ccledata$copies)) = c("gene", "CCLE_ID")
names(dimnames(ccledata$eset)) = c("gene", "CCLE_ID")
cd = ccledata$clines %>%
    select(CCLE_ID, Name, Site_Primary) %>%
    left_join(reshape2::melt(ccledata$copies[top,], value.name="copies")) %>%
    left_join(reshape2::melt(ccledata$eset[top,], value.name="expr")) %>%
    mutate(expr = expr * copies/2) %>% # undo normmatrix normalization
    group_by(gene) %>%
        filter(expr < quantile(expr, 0.95)) %>%
        filter(copies < quantile(copies, 0.95)) %>%
    ungroup()
abl = cd %>%
    group_by(gene) %>%
    summarize(mean = median(expr[abs(copies-2) < 1.8]))
pccle = ggplot(cd, aes(x=copies, y=expr)) +
    annotate("rect", xmin=1.8, xmax=2.2, ymin=-Inf, ymax=Inf, alpha=0.2, fill="yellow") +
    geom_vline(xintercept=2, color="grey") +
    geom_vline(xintercept=c(1.8,2.2), color="grey", linetype="dotted") +
    geom_abline(data=abl, aes(intercept=0, slope=mean/2), color="red", linetype="dashed") +
    geom_point(alpha=0.3) +
    geom_smooth(method="lm") +
    facet_wrap(~ gene, scales="free") +
    labs(title = "CCLE compensation (red: expected, blue: observed); 95th% shown (expr/copies); yellow=euploid",
         y = "normalized read count")

pdf(args$plotfile, 14, 12)
cowplot::plot_grid(plotlist=overview)
print(porf)
print(pccle)
dev.off()
