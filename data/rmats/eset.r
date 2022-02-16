library(dplyr)
library(ggplot2)
sys = import('sys')
plt = import('plot')

#' Plot read count distributions
#'
#' @param reads  Rsubread::featureCounts object
#' @return       ggplot2 object with bars for total, mapped, assigned reads
read_count_plot = function(reads) {
    stats = reads$stat %>%
        tidyr::gather("bam", "value", -Status) %>%
        group_by(bam) %>%
        summarize(total = sum(value),
                  mapped = sum(value[!grepl("Unmapped", Status)]),
                  assigned = value[Status == "Assigned"]) %>%
        mutate(bam = sub("(\\.sorted)?\\.bam$", "", bam))

    long = stats %>%
        tidyr::gather("type", "value", -bam) %>%
        mutate(type = factor(type, levels=c("total", "mapped", "assigned")))

    ggplot(long, aes(x=forcats::fct_reorder(bam, -value), y=value)) +
        geom_col(aes(fill=type), width=2, position=position_dodge(width=0.3)) +
        scale_fill_manual(values=c(total="#cecece", mapped="#00aedb55", assigned="#d11141ce")) +
        theme_minimal() +
        coord_cartesian(expand=FALSE) +
        theme(axis.text.x = element_text(angle=45, hjust=1),
              axis.title.x = element_blank())
}

#' Plot a PCA of samples
#'
#' @param eset  DESeq2 transform or expression set object
#' @return      ggplot2 object
plot_pca = function(eset, mapping=aes(x=PC1, y=PC2)) {
    shapes = c("8h"=21, "24h"=22, "72h"=23)
    sizes = c(Luc=2, RBM14=4)
    cols = RColorBrewer::brewer.pal(length(levels(eset$cline)), "Set1") %>% setNames(levels(eset$cline))

    plt$pca(eset, mapping) +
        geom_point(aes(fill=cline, shape=time, size=cond), alpha=0.7) +
        scale_shape_manual(values=shapes, guide=guide_legend(override.aes=list(size=3))) +
        scale_size_manual(values=sizes) +
        scale_fill_manual(values=cols, guide=guide_legend(override.aes=list(shape=21, size=3))) +
        theme_minimal() +
        ggrepel::geom_text_repel(aes(label=paste(sprintf("%s:%s:%s:%i", cline, cond, time, rep))))
}

plot_rbm14 = function(eset) {
    eset = DESeq2::estimateSizeFactors(eset)
    dset = cbind(as.data.frame(SummarizedExperiment::colData(eset)),
                 RBM14=DESeq2::counts(eset, normalized=TRUE)["RBM14",]) %>%
        mutate(label = paste(cline, cond, time, rep, sep=":"),
               group = paste(cline, cond, rep, sep=":"))

    ggplot(dset, aes(x=time, y=RBM14, group=group, color=cond)) +
        geom_line(aes(linetype=rep), size=1.5, alpha=0.7) +
        facet_wrap(~ cline, ncol=1, scales="free") +
        theme_minimal() +
        ylim(c(0, NA)) +
        scale_color_manual(values=c(Luc="grey", RBM14="red")) +
        labs(y = "RBM14 (normalized read count)")
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'eset.rds'),
        opt('p', 'plotfile', 'rds', 'eset.pdf')
    )

    bams = list.files("seqdata", "\\.bam$", recursive=TRUE, full.names=TRUE)
    reads = Rsubread::featureCounts(bams, annot.ext="seqdata/hg38.refseq.gtf",
                    isGTFAnnotationFile=TRUE, isPairedEnd=TRUE, nthreads=10)

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
    vst = DESeq2::varianceStabilizingTransformation(eset)

    pdf(args$plotfile, 10, 8)
    print(read_count_plot(reads))
    print(plot_pca(vst))
    print(plot_pca(vst, aes(x=PC3, y=PC4)))
    print(plot_pca(vst[,vst$time == "8h"], aes(x=PC3, y=PC4)) + ggtitle("8h"))
    print(plot_pca(vst[,vst$time == "24h"], aes(x=PC3, y=PC4)) + ggtitle("24h"))
    print(plot_pca(vst[,vst$cline == "H838"]))
    print(plot_pca(vst[,vst$cline == "H1650"]))
    print(plot_pca(vst[,vst$cline == "HCC70"]))
    print(plot_rbm14(eset))
    dev.off()

    saveRDS(eset, file=args$outfile)
})
