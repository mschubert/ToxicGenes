library(dplyr)
library(patchwork)
sys = import('sys')
plt = import('plot')

venn = function(title, df, comp=-0.5, fdr=c(0.05,1e-3), r2=0.02) {
    df2 = df %>%
        filter(estimate < comp,
               (dset %in% c("orf","ccle") & adj.p < fdr[1]) |
                    (dset == "tcga" & adj.p < fdr[2]),
               dset=="orf" | rsq > r2) %>%
        select(name, dset) %>% distinct()
    df3 = unstack(df2)
    df3 = df3[sapply(df3, length) > 0]
    p1 = plt$try(plt$venn(df3)) +
        ggtitle(sprintf("%s >= %i%% comp, %i/%.2g%% FDR, %i%% R^2", title,
                        round(-comp*100), round(fdr[1]*100), fdr[2]*100, round(r2*100)))
    df4 = df2 %>%
        group_by(name) %>%
        mutate(n = n()) %>%
        ungroup() %>%
        arrange(-n) %>%
        filter(n >= 2)
    if (length(unique(df4$name)) > 200)
        df4 = df4 %>% filter(n >= 3)
    df4 = df4 %>%
        mutate(name = factor(name, levels=rev(sort(unique(name)))))
    tsize = min(6, round(450 / length(levels(df4$name))))
    p2 = ggplot(df4, aes(x=factor(dset, levels=c("orf", "ccle", "tcga")), y=name)) +
        geom_tile(aes(fill=factor(n)), color="white", size=0.5, show.legend=FALSE) +
        coord_fixed() +
        scale_fill_manual(values=c("black", "red")) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.x = element_text(size=tsize, angle=90, hjust=1, vjust=0.5),
              axis.text.y = element_text(size=tsize))
    p1 + p2 + plot_layout(widths=c(5,1))
}

sys$run({
    args = sys$cmd$parse(
        opt('d', 'dset', 'merge', '../merge/pan.rds'),
        opt('f', 'fit', 'rank|rlm{,2,3}', 'rlm3'),
        opt('p', 'plotfile', 'pdf', 'pan_rlm3.pdf'))

    dset = readRDS(args$dset) %>%
        filter(dset == "orf" | fit == args$fit)

    if (args$fit == "rlm2")
        dset$rsq = 1 # not implemented in rlm2 model, do not filter on it

    plots = list(
        venn("amp+del", dset %>% filter(cna %in% c("oe", "all"))),
        venn("amp+del", dset %>% filter(cna %in% c("oe", "all")), fdr=c(0.05,0.01), r2=0.1),
        venn("amp", dset %>% filter(cna %in% c("oe", "amp"))),
        venn("amp", dset %>% filter(cna %in% c("oe", "amp")), fdr=c(0.05,0.01), r2=0.1),
        venn("del", dset %>% filter(cna %in% c("del")))
    )
    pdf(args$plotfile)
    for (p in plots)
        print(p)
    dev.off()
})
