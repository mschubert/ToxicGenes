io = import('io')
sys = import('sys')
idmap = import('process/idmap')
seq = import('seq')

args = sys$cmd$parse(
    opt('a', 'annot', 'txt', 'Cell_lines_annotations_20181226.txt'),
    opt('r', 'rnaseq', 'gct', 'CCLE_RNAseq_genes_counts_20180929.gct.gz'),
    opt('c', 'copies', 'txt', 'CCLE_copynumber_byGene_2013-12-03.txt.gz'),
    opt('o', 'outfile', 'rds', 'dset.rds'))

idx = readr::read_tsv(args$annot)

expr = io$read_gct(args$rnaseq)@mat
expr = expr[rowMeans(expr) >= 10,] # small bias for library size?
rownames(expr) = idmap$gene(sub("\\.[0-9]+", "", rownames(expr)), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)),]

cinfo = readr::read_tsv(args$copies) # assuming log2(x/2)
copies = 2^data.matrix(cinfo[,6:ncol(cinfo)]) + 1
rownames(copies) = cinfo$SYMBOL

chrs = seq$coords$gene(chromosomes=1:22)

narray::intersect(idx$CCLE_ID, expr, copies, along=2)
narray::intersect(chrs$external_gene_name, expr, copies, along=1)

saveRDS(list(idx=idx, expr=expr, copies=copies), file=args$outfile)
