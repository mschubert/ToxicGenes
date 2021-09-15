library(dplyr)
sys = import('sys')
gset = import('data/genesets')
plt = import('plot')

test_set = function(dset, sets, set, cna) {
    cnv = dset[[cna]]
    in_set = cnv$gene %in% sets[[set]]
    lm(estimate ~ in_set, data=cnv) %>% # statistic if rlm, 'estimate' if stan-nb
        broom::tidy() %>%
        filter(term == "in_setTRUE") %>%
        select(-term) %>%
        mutate(size = sum(in_set),
               genes = ifelse(sum(in_set) <= 10, paste(sets[[set]], collapse=","), NA_character_))
}

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('i', 'infile', 'xlsx', 'pan/stan-nb.xlsx'),
    opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
    opt('p', 'plotfile', 'pdf', 'CH.HALLMARK.pdf')
)

cnas = yaml::read_yaml(args$config)$cna

dset = readxl::excel_sheets(args$infile) %>%
    sapply(function(s) readxl::read_xlsx(args$infile, sheet=s), simplify=FALSE)
names(dset) = "amp" #FIXME:

sets = readRDS(args$setfile) %>%
    gset$filter(min=1, valid=dset$all$name)

result = expand.grid(cna=cnas, set=names(sets), stringsAsFactors=FALSE) %>%
    mutate(result = clustermq::Q_rows(., test_set, const=list(dset=dset, sets=sets), n_jobs=0)) %>%
    tidyr::unnest(cols="result") %>%
    group_by(cna) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"),
           set = ifelse(is.na(genes), set, sprintf("%s\n%s", set, genes))) %>%
    arrange(adj.p, p.value) %>%
    ungroup()

do_plot = . %>%
    mutate(label = set) %>%
    plt$volcano(base.size=0.2, text.size=2.5, label_top=20, repel=TRUE,
                pos_label_bias=0.2)

plots = result %>%
    split(.$cna) %>%
    lapply(do_plot)

pdf(args$plotfile, 10, 8)
for (i in seq_along(plots))
    print(plots[[i]] + ggtitle(names(plots)[i]))
dev.off()
