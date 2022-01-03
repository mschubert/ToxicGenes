library(dplyr)
library(ggplot2)
sys = import('sys')

read_count_plot = function(eset) {
    df = as.data.frame(SummarizedExperiment::colData(eset)) %>%
        mutate(name = sub("(\\.sorted)?\\.bam", "", colnames(eset)),
               counts = colSums(SummarizedExperiment::assay(eset)))
    ggplot(df, aes(x=forcats::fct_reorder(name, -counts), y=counts)) +
        geom_col() +
        theme(axis.text.x = element_text(angle=45, hjust=1))
}

plot_pca = function(eset) {
    vst = DESeq2::varianceStabilizingTransformation(eset)

    pcadata = DESeq2::plotPCA(vst, intgroup=c("cline", "cond", "time", "rep"), returnData=TRUE)
    pcavar = round(100 * attr(pcadata, "percentVar"))
    shapes = c("8h"=21, "24h"=22, "72h"=23)
    sizes = c("Luc"=2, "RBM14"=4)

    ggplot(pcadata, aes(x=PC1, y=PC2)) +
        geom_point(aes(fill=cline, shape=time, size=cond), alpha=0.7) +
        scale_shape_manual(values=shapes) +
        scale_size_manual(values=sizes) +
        xlab(paste0("PC1: ", pcavar[1], "% variance")) +
        ylab(paste0("PC2: ", pcavar[2], "% variance")) +
        ggrepel::geom_text_repel(aes(label=group))
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'eset.rds'),
        opt('p', 'plotfile', 'rds', 'eset.pdf')
    )

    bams = list.files("seqdata", "\\.bam$", recursive=TRUE, full.names=TRUE)
    reads = Rsubread::featureCounts(bams, annot.ext="seqdata/ref_annot.gtf",
                                    isGTFAnnotationFile=TRUE, isPairedEnd=TRUE)

    meta = basename(bams) %>%
        sub(".sorted.bam", "", ., fixed=TRUE) %>%
        strsplit("_") %>%
        do.call(rbind, .) %>%
        as.data.frame() %>%
        transmute(cline = factor(V1),
                  cond = sub("luc", "Luc", V2) %>% factor() %>% relevel("Luc"),
                  time = factor(V3) %>% relevel("8h"),
                  rep = factor(V4)) %>%
        as_tibble()

    eset = DESeq2::DESeqDataSetFromMatrix(reads$counts, meta, ~1)

    pdf(args$plotfile)
    print(read_count_plot(eset))
    print(plot_pca(eset))
    dev.off()

    saveRDS(eset, file=args$outfile)
})
