library(dplyr)
sys = import('sys')
gset = import('genesets')
plt = import('plot')

args = sys$cmd$parse(
    opt('i', 'infile', 'xlsx', 'pan/betareg.xlsx'),
    opt('s', 'setfile', 'rds', '../data/genesets/MSigDB_Hallmark_2020.rds'),
    opt('p', 'plotfile', 'pdf', 'pan/betareg/MSigDB_Hallmark_2020.pdf')
)

dset = readxl::read_xlsx(args$infile)
sets = readRDS(args$setfile) %>%
    gset$filter(min=1, valid=dset$gene)

res = gset$test(dset, sets)

pdf(args$plotfile, 10, 8)
print(plt$volcano(res, base.size=0.2, text.size=2.5, label_top=20, pos_label_bias=0.2))
dev.off()
