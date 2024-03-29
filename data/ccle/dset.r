io = import('io')
sys = import('sys')
idmap = import('process/idmap')
seq = import('seq')

args = sys$cmd$parse(
    opt('a', 'annot', 'txt', 'Cell_lines_annotations_20181226.txt'),
    opt('r', 'rnaseq', 'gct', 'CCLE_RNAseq_genes_counts_20180929.gct.gz'),
    opt('c', 'copies', 'txt', 'CCLE_copynumber_byGene_2013-12-03.txt.gz'),
    opt('s', 'mut', 'csv', 'CCLE_mutations.csv'),
    opt('m', 'meth', 'txt', 'CCLE_RRBS_TSS_1kb_20180614.txt'),
    opt('o', 'outfile', 'rds', 'dset.rds'))

genes = seq$coords$gene(c("ensembl_gene_id", "external_gene_name"), chromosomes=1:22)
clines = readr::read_tsv(args$annot) %>%
    mutate(tcga_code = sub("/", "", tcga_code),
           tcga_code = ifelse(tcga_code == "UNABLE TO CLASSIFY", NA, tcga_code))

###
### Gene expression
###
expr2 = io$read_gct(args$rnaseq)@mat
hgnc = idmap$gene(sub("\\.[0-9]+", "", rownames(expr2)), to="hgnc_symbol")
uhgnc = unique(na.omit(hgnc))
expr = matrix(NA_integer_, nrow=length(uhgnc), ncol=ncol(expr2),
              dimnames = list(uhgnc, colnames(expr2)))
for (n in rownames(expr))
    expr[n,] = colSums(expr2[which(n == hgnc),,drop=FALSE])
stopifnot(sum(is.na(expr)) == 0)
#expr = narray::map(expr, along=1, sum, subsets=rownames(expr)) # oom
expr = expr[rowMeans(expr) >= 10,] # small library size bias

###
### Copy number segments -> genes
###
cinfo = readr::read_tsv(args$copies) # assuming y = log2(x)-1 -> x = 2^(y+1)
copies = 2 ^ (data.matrix(cinfo[,6:ncol(cinfo)]) + 1)
e2hgnc = idmap$gene(as.character(cinfo$EGID), from="entrezgene_id", to="external_gene_name")
pick = function(x1, x2) { # 16k -> 17k genes recovered where we have expr
    if (x1 %in% rownames(expr)) { x1
    } else if (x2 %in% rownames(expr)) { x2
    } else if (x1 %in% genes$external_gene_name) { x1
    } else if (x2 %in% genes$external_gene_name) { x2
    } else { x1 }
}
rownames(copies) = mapply(pick, x1=cinfo$SYMBOL, x2=e2hgnc)

###
### Mutations
###
mut = readr::read_csv(args$mut)  %>%
    inner_join(clines %>% select(cline=Name, DepMap_ID=depMapID)) %>%
    transmute(cline = cline,
              gene = Hugo_Symbol,
              type = Variant_Classification)

###
### Methylation
###
mdata = readr::read_tsv(args$meth) 
meth_all = data.matrix(mdata[,8:ncol(mdata)])
rownames(meth_all) = mdata$gene

###
### Merge, add LOC from ORF screen
###
narray::intersect(clines$CCLE_ID, expr, copies, along=2)
narray::intersect(genes$external_gene_name, expr, copies, along=1)
meth = matrix(NA, ncol=ncol(expr), nrow=nrow(expr), dimnames=dimnames(expr))
common_row = intersect(rownames(meth), rownames(meth_all))
common_col = intersect(colnames(meth), colnames(meth_all))
meth[common_row, common_col] = meth_all[common_row, common_col]

eset_raw = DESeq2::DESeqDataSetFromMatrix(expr, clines, ~1) %>%
    DESeq2::estimateSizeFactors(normMatrix=copies)
eset = DESeq2::counts(eset_raw, normalized=TRUE)

#loc = readr::read_tsv("LOC254896_GE_values.txt") %>%
#    left_join(readr::read_csv("./CCLE_DepMapID_name_mapping.csv") %>%
#                 rename(DepMap_ID = broad_id))
#eset = rbind(eset, LOC254896=loc$LOC254896[match(colnames(eset), loc$ccle_name)])
#copies = rbind(copies, LOC254896=copies["TNFRSF10C",])

saveRDS(list(clines=clines, copies=copies, eset_raw=eset_raw, eset=eset, mut=mut, meth=meth), file=args$outfile)
