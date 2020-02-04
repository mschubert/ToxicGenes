library(dplyr)
sys = import('sys')
gset = import('data/genesets')
plt = import('plot')

test_set = function(dset, sets, set, cna) {
    cnv = dset[[cna]]
    in_set = cnv$name %in% sets[[set]]
    lm(statistic ~ in_set, data=cnv) %>%
        broom::tidy() %>%
        filter(term == "in_setTRUE") %>%
        select(-term) %>%
        mutate(size = sum(in_set))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'xlsx', 'pan/rlm3_pur.xlsx'),
    opt('s', 'setfile', 'rds', '../data/genesets/GO_Biological_Process_2018.rds'),
    opt('p', 'plotfile', 'pdf', 'GO_Biological_Process_2018.pdf'))

dset = readxl::excel_sheets(args$infile) %>%
    sapply(function(s) readxl::read_xlsx(args$infile, sheet=s), simplify=FALSE)

sets = readRDS(args$setfile) %>%
    gset$filter(min=4, valid=dset$all$name)

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
    plt$volcano(base.size=0.2, text.size=2.5, label_top=30, repel=TRUE)

plots = result %>%
    split(.$cna) %>%
    lapply(do_plot)

pdf(args$plotfile)
for (i in seq_along(plots))
    print(plots[[i]] + ggtitle(names(plots)[i]))
dev.off()
