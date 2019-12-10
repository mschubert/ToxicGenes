library(dplyr)
sys = import('sys')
gset = import('data/genesets')
plt = import('plot')

test_set = function(dset, sets, set, cna) {
    dset = dset[dset$cna %in% c(cna, "oe"),]
    in_set = dset$name %in% sets[[set]]
    lm(statistic ~ in_set, data=dset) %>%
        broom::tidy() %>%
        filter(term == "in_setTRUE") %>%
        select(-term) %>%
        mutate(size = sum(in_set))
}

args = sys$cmd$parse(
    opt('d', 'dset', 'merged rds', '../merge/pan.rds'),
    opt('f', 'fit', 'rank|rlm{,2,3}', 'rlm3'),
    opt('s', 'sets', 'rds', '../data/genesets/GO_Biological_Process_2018.rds'),
    opt('p', 'plotfile', 'pdf', 'GO_Biological_Process_2018.pdf'))

dset = readRDS(args$dset) %>%
    filter(fit %in% c("lm", args$fit),
           dset != "tcga" | adj == "pur")

sets = readRDS(args$sets) %>%
    gset$filter(min=4, valid=dset$name)

result = expand.grid(cna=c("all", "amp", "del"), set=names(sets), stringsAsFactors=FALSE) %>%
    mutate(result = clustermq::Q_rows(., test_set, const=list(dset=dset, sets=sets), n_jobs=0)) %>%
    tidyr::unnest(cols="result") %>%
    group_by(cna) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value) %>%
    ungroup()

do_plot = . %>%
    mutate(label = set) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
    plt$volcano(text.size=2.5, label_top=30, repel=TRUE)

plots = result %>%
    split(.$cna) %>%
    lapply(do_plot)

pdf(args$plotfile)
for (i in seq_along(plots))
    print(plots[[i]] + ggtitle(names(plots)[i]))
dev.off()
