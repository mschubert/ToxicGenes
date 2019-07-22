io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'dset', 'merge_genes.rds'),
    opt('s', 'select', 'genes|gene set', 'genes'),
    opt('n', 'num', 'number of top hits', '10'),
    opt('o', 'outfile', 'yaml', 'top_genes.yaml'))

dset = readRDS(args$infile)

top = dset %>%
    group_by(dset, fit, adj) %>%
    mutate(score = rank(-statistic) / length(statistic)) %>%
    group_by(name, dset, fit) %>%
    summarize(score = mean(score, na.rm=TRUE)) %>%
    group_by(name, dset) %>%
    summarize(score = mean(score, na.rm=TRUE)) %>%
    group_by(name) %>%
    filter(all(score > 0)) %>%
    summarize(n_dset = length(name),
              score = sum(score, na.rm=TRUE) / n_dset^0) %>% # ^0 is sum (CDKN1A), ^1 is mean (LOCxxx)
    arrange(-score)

obj = list(head(top$name, as.integer(args$num)))
names(obj) = args$select
yaml::write_yaml(obj, file=args$outfile)
