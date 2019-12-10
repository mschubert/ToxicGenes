library(dplyr)
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'pan merge', '../merge/pan.rds'),
    opt('p', 'plotfile', 'pdf', 'rank_overview.pdf'))

dset = readRDS(args$dset)
