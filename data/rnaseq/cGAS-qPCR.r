library(dplyr)
library(ggplot2)
library(DESeq2)

eset = readRDS("eset.rds") |>
    DESeq2::estimateSizeFactors()
genes = c("IL6", "CXCL8", "CXCL10", "CCL5")

dset = t(counts(eset, normalized=TRUE)[genes,]) |>
    cbind(as.data.frame(colData(eset))) |>
    tidyr::pivot_longer(all_of(genes))

pdf("cGAS-qPCR.pdf")
ggplot(dset, aes(x=cline, fill=cond, y=value)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(aes(shape=time), size=2, alpha=0.8) +
    scale_shape_manual(values=c(21,22,23)) +
    facet_wrap(~name, scales="free") +
    scale_y_continuous(trans="log1p")
dev.off()
