library(dplyr)
library(magrittr)
library(DESeq2)
sys = import('sys')
plt = import('plot')
gset = import('genesets')

plot_pca = function(vst) {
    pcadata = DESeq2::plotPCA(vst, intgroup=colnames(colData(vst)), returnData=TRUE)
    pcavar = round(100 * attr(pcadata, "percentVar"))

    ggplot(pcadata, aes(x=PC1, y=PC2)) +
        geom_point(aes(color=treatment, shape=time), alpha=0.9, size=3) +
        ggrepel::geom_text_repel(aes(label=sample)) +
        xlab(paste0("PC1: ", pcavar[1], "% variance")) +
        ylab(paste0("PC2: ", pcavar[2], "% variance"))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'csv', 'STAR_Gene_Counts.csv'),
    opt('o', 'outfile', 'rds', 'dset_both.rds'),
    opt('p', 'plotfile', 'pdf', 'dset_both.pdf')
)

dset = readr::read_csv(args$infile)
emat = data.matrix(dset[2:13])
rownames(emat) = dset$Gene_ID
meta = strsplit(colnames(emat), "_") %>% do.call(rbind, .) %>%
    as.data.frame() %>% as_tibble() %>%
    dplyr::transmute(cell_line = V1,
                     treatment = factor(V2) %>% relevel("luc"),
                     time = factor(V3) %>% relevel("8h"),
                     rep = factor(V4),
                     sample = paste(treatment, time, rep))

eset = DESeq2::DESeqDataSetFromMatrix(emat, meta, ~1)
vs = DESeq2::varianceStabilizingTransformation(eset)
colData(eset)$RBM14 = assay(vs)["RBM14",]


pdf(args$plotfile)

plot_pca(vs)

ggplot(as.data.frame(colData(eset)), aes(x=time, y=RBM14)) +
    geom_point(aes(color=treatment), size=5) +
    facet_wrap(~ treatment)

dev.off()

saveRDS(eset, file=args$outfile)
