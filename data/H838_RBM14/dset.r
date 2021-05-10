library(dplyr)
library(magrittr)
library(DESeq2)
sys = import('sys')
plt = import('plot')
gset = import('genesets')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../../config.yaml'),
#    opt('i', 'infile', 'txt', './210331_VR8768_viper/STAR/H838_luc_1/H838_luc_1.counts.tab'),
    opt('o', 'outfile', 'rds', 'dset_both.rds'),
    opt('p', 'plotfile', 'pdf', 'dset_both.pdf')
)

cfg = yaml::read_yaml(args$config)

smps = c("H838_luc_1", "H838_luc_2", "H838_RBM14_1", "H838_RBM14_2")
stypes = sub(".*(luc|RBM14).*", "\\1", smps)
reads = sprintf("./210331_VR8768_viper/STAR/%s/%s.counts.tab", smps, smps) %>%
    lapply(readr::read_tsv, skip=4, col_names=FALSE) %>%
    lapply(. %$% setNames(X2, X1)) %>%
    setNames(smps) %>%
    narray::stack(along=2)
reads = reads[rowSums(reads) != 0,]

meta = data.frame(condition = factor(stypes, levels=c("luc", "RBM14")))
eset = DESeq2::DESeqDataSetFromMatrix(reads, meta, ~condition)
res = DESeq2::DESeq(eset) %>%
    DESeq2::results(name = "condition_RBM14_vs_luc") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_name") %>%
    as_tibble() %>%
    arrange(padj, pvalue)

res$size = 1 # fixme @ volcano util

#fixme: gset$get_human better error on invalid set name
sets = sapply(cfg$genesets, function(s) readRDS(paste0("../genesets/", s, ".rds")), simplify=FALSE)
scores = lapply(sets, gset$test, genes=res)

pdf(args$plotfile, 10, 8)
print(plt$volcano(res) + ggtitle("genes"))
for (i in seq_along(scores))
    print(plt$volcano(scores[[i]]) + ggtitle(names(scores)[i]))
dev.off()

saveRDS(eset, file=args$outfile)
