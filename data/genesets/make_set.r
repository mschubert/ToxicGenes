sys = import('sys')
idmap = import('process/idmap')
enr = import('tools/enrichr')
msdb = import('tools/msigdb')

args = sys$cmd$parse(
    opt('g', 'geneset', 'Identifier of the gene set', 'KEA_2015'),
    opt('o', 'outfile', 'save to rds', 'KEA_2015.RData'))

if (args$geneset %in% enr$dbs()$name) {
    sets = enr$genes(args$geneset)
} else if (args$geneset %in% msdb$dbs()) {
    sets = msdb$genes(args$geneset)
} else
    stop("invalid gene set: ", args$geneset)

saveRDS(sets, file=args$outfile)
