library(dplyr)
library(ggplot2)
sys = import('sys')

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'eset_exon.rds')
    )

    bams = list.files("seqdata", "\\.bam$", recursive=TRUE, full.names=TRUE)
    reads = Rsubread::featureCounts(bams, annot.ext="seqdata/hg38.refseq.gtf",
                    isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, nthreads=10,
                    allowMultiOverlap=TRUE, juncCounts=TRUE, useMetaFeatures=FALSE)
    rownames(reads$counts) = make.unique(rownames(reads$counts))

    meta = sub(".sorted.bam", "", basename(bams), fixed=TRUE) %>%
        strsplit("_") %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        transmute(cline = factor(V1),
                  cond = sub("luc", "Luc", V2) %>% factor() %>% relevel("Luc"),
                  time = factor(V3) %>% relevel("8h"),
                  rep = factor(V4)) %>%
        as_tibble()

    eset = DESeq2::DESeqDataSetFromMatrix(reads$counts, meta, ~1)
    eset = eset[rowSums(SummarizedExperiment::assay(eset)) != 0,]

    saveRDS(eset, file=args$outfile)
})
