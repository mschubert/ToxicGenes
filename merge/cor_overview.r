library(dplyr)
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'pan merge', '../merge/pan.rds'),
    opt('p', 'plotfile', 'pdf', 'cor_overview.pdf'))

dset = readRDS(args$dset)
