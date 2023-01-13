library(dplyr)
library(SummarizedExperiment)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
idmap = import('process/idmap')
deseq = import('process/deseq')
dori = import('../../../dorine/data/rnaseq_adaptation/counts')

# merge count files present here
reads = list.files("aligned", "ReadsPerGene\\.out\\.tab", full.names=TRUE) %>%
    dori$read_files() %>%
    mutate(sample = sub(".*_([0-9A-Z]+)$", "\\1", sample))
stats = dori$read_stats(reads)

# do standard qc plots
sel = stats %>%
    group_by(count_type) %>%
        summarize(part_mapped = sum(reads[feature == "N_mapped"]) / sum(reads)) %>%
    slice_max(part_mapped) %>%
    pull(count_type)
levels(stats$count_type)[levels(stats$count_type) == sel] = sprintf("%s (selected)", sel)
dori$read_barplot(stats)

cmat = dori$read_matrix(reads, sel)
meta = data.frame(sample=colnames(cmat))
eset = DESeq2::DESeqDataSetFromMatrix(cmat, meta, ~1)
plt$pca(eset) +
    geom_point() +
    geom_text(aes(label=sample))

# add dorine's RPE1s to combined DESeq data set
nki = readRDS("../../../dorine/data/rnaseq_adaptation/eset.rds")
ref = nki[,nki$sample == "p53kd"]
common = intersect(rownames(eset), rownames(ref))
colData(ref) = colData(ref)["sample"]
both = cbind(eset[common,], ref[common,])
plt$pca(both) +
    geom_point() +
    geom_text(aes(label=sample))

# compute DE over the ref sample
both$sample = relevel(factor(both$sample), "RPE1SS48")
res = deseq$genes(both, ~sample, extract="^sample_RPE") %>%
    mutate(term = sub("sample_RPE1(.+)_vs_.*", "\\1", term))

