io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'dset', 'merge/pan/genes.rds'),
    opt('s', 'select', 'genes|gene set', 'genes'),
    opt('n', 'num', 'number of top hits', '12'),
    opt('o', 'outfile', 'yaml', 'pan/top.stat_genes.yaml'))

dset = readRDS(args$infile)

top = dset %>%
    filter(adj %in% c("none", "puradj"),
           fit %in% c("lm", "rank", "rlm")) %>%
    group_by(dset, fit, adj) %>%
    mutate(score = (1-adj.p) * rank(-statistic) / length(statistic)) %>%
    group_by(name, dset, fit) %>% # select most significant cna (amp, del, all)
    top_n(1, -score) %>%
    group_by(name, dset) %>% # mean by fit (rlm, rank)
    summarize(score = mean(score, na.rm=TRUE)) %>%
    group_by(name) %>% # summarize by dset (orf, ccle, tcga)
    summarize(n_dset = length(name),
              score = sum(score^0.5, na.rm=TRUE)) %>% # root prioritizes consistency across dsets
    arrange(-score)

obj = list(head(top$name, as.integer(args$num)))
names(obj) = args$select
obj$methods = c("lm", "rank", "rlm")
yaml::write_yaml(obj, file=args$outfile)
