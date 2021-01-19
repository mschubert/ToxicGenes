library(dplyr)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'ccle/dset.rds'),
    opt('o', 'outfile', 'xlsx', 'df_ccle.rds')
)

dset = readRDS(args$infile)

eset = DESeq2::estimateSizeFactors(dset$eset_raw)
copies = dset$copies

cnts = DESeq2::counts(eset)
names(dimnames(cnts)) = names(dimnames(copies)) = c("gene", "cell_line")
sfs = data.frame(cell_line=colnames(eset), sf=DESeq2::sizeFactors(eset))
cv = data.frame(cell_line=colnames(eset), covar=dset$clines$tcga_code)
cnts = reshape2::melt(cnts, value.name="expr")
copies2 = reshape2::melt(copies, value.name="copies")

df = inner_join(cnts, copies2) %>%
    inner_join(sfs) %>%
    inner_join(cv) %>%
    na.omit() %>%
    mutate(eup_equiv = (copies - 2) / 2)

saveRDS(df, file=args$outfile)
