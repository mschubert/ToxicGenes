io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'dset', 'merge/pan/genes.rds'),
    opt('s', 'select', 'genes|gene set', 'genes'),
    opt('n', 'num', 'number of top hits', '12'),
    opt('o', 'outfile', 'yaml', 'pan/top.stat_genes.yaml'))

n = as.integer(args$num)
dset = readRDS(args$infile) %>%
    filter(adj %in% c("none", "puradj"),
           fit %in% c("lm", "rank", "rlm")) %>%
    mutate(score = (1-adj.p) * rank(-statistic) / length(statistic))

select_top = . %>%
    group_by(name, dset) %>% # mean by fit (rlm, rank)
    summarize(score = mean(score, na.rm=TRUE)) %>%
    group_by(name) %>% # summarize by dset (orf, ccle, tcga)
    summarize(n_dset = length(name),
              score = sum(score^0.5, na.rm=TRUE)) %>% # root prioritizes consistency across dsets
    arrange(-score)

top_all = select_top(dset %>% filter(cna %in% c("oe", "all")))
top_amp = select_top(dset %>% filter(cna %in% c("oe", "amp")))
top_del = select_top(dset %>% filter(cna %in% c("oe", "del")))

obj = list(list(
    all = top_all$name[seq_len(n)],
    all2 = top_all$name[seq_len(n)+n],
    amp = top_amp$name[seq_len(n)],
    del = top_del$name[seq_len(n)]
))
names(obj) = args$select
obj$methods = c("lm", "rank", "rlm")
yaml::write_yaml(obj, file=args$outfile)
