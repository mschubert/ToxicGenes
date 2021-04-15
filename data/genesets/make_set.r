library(dplyr)
sys = import('sys')
gset = import('genesets')

args = sys$cmd$parse(
    opt('g', 'geneset', 'Identifier of the gene set', 'MSigDB_Hallmark_2020'),
    opt('o', 'outfile', 'rds', 'MSigDB_Hallmark_2020.rds')
)

sets = gset$get_human(args$geneset)
saveRDS(sets, file=args$outfile)
