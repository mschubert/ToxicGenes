library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'plot_cnvs.rds'),
    opt('m', 'meta', 'xlsx', 'N4plus_data_Michael_20230414.xlsx'),
    opt('p', 'plotfile', 'pdf', 'surv.pdf')
)

meta = readxl::read_xlsx(args$meta) %>%
    mutate_all(function(x) ifelse(x == 9999, NA, x)) %>%
    mutate(PtID = sprintf("%03i", PtID))

dset = readRDS(args$infile)

