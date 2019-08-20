io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'dset', 'merge/pan/genes.rds'),
    opt('s', 'select', 'genes|gene set', 'genes'),
    opt('n', 'num', 'number of top hits', '12'),
    opt('o', 'outfile', 'yaml', 'pan/top.stat_genes.yaml'))

n = as.integer(args$num)
meths = c("lm", "rlm", "rlm2")

dset = readRDS(args$infile) %>%
    filter(adj %in% c("none", "puradj"),
           fit %in% c("lm", "rlm2")) %>%
    mutate(estimate = ifelse(dset == "orf", 2^estimate - 1, estimate),
           estimate = pmax(estimate, -1)) %>%
    mutate(score = (1-adj.p) * (rank(-statistic) / length(statistic)) * (-estimate))

select_top = . %>%
    group_by(name, dset) %>% # mean by fit (rlm, rank)
    summarize(score = mean(score, na.rm=TRUE)) %>%
    group_by(name) %>% # summarize by dset (orf, ccle, tcga)
    summarize(n_dset = length(name),
              score = sum(score^0.5, na.rm=TRUE)) %>% # root prioritizes consistency across dsets
    arrange(-score)

top_all = select_top(dset %>% filter(cna %in% c("oe", "all"))) %>% pull(name) %>% head(n*2)
top_amp = select_top(dset %>% filter(cna %in% c("oe", "amp"), ! name %in% top_all)) %>% pull(name)
top_del = select_top(dset %>% filter(cna %in% c("oe", "del"), ! name %in% top_all)) %>% pull(name)

obj = list(list(
    "top 12 (amp and del)" = top_all[seq_len(n)],
    "top 24 (amp and del)" = top_all[seq_len(n)+n],
    "amp top 12" = top_amp[seq_len(n)],
    "del top 12" = top_del[seq_len(n)]
))
names(obj) = args$select
obj$methods = meths
yaml::write_yaml(obj, file=args$outfile)
