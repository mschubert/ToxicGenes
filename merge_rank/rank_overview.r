library(dplyr)
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'pan merge', '../merge/pan.rds'),
    opt('p', 'plotfile', 'pdf', 'rank_overview.pdf'))

dset = readRDS(args$dset) %>%
    mutate(dset = paste(dset, cna, fit, adj)) %>%
    select(name, dset, statistic)

mat = narray::construct(statistic ~ name + dset, data=dset)
cmat = cor(na.omit(mat))
