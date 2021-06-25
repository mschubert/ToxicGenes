library(dplyr)
library(GenomicRanges)
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')

plot_ccle_meth = function(cd, et=0.15) {
    mut_shapes = c(
        Missense_Mutation = 23,
        Nonsense_Mutation = 25,
        Frame_Shift_Del = 25,
        Frame_Shift_Ins = 25,
        Splice_Site = 24,
        In_Frame_Del = 24,
        Silent = 22
    )

    ggplot(cd, aes(x=copies, y=expr)) +
        annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, alpha=0.2, fill="yellow") +
        geom_vline(xintercept=2, color="grey") +
        geom_vline(xintercept=c(2-et,2+et,1+et,3-et), color="grey", linetype="dotted") +
        geom_point(aes(shape=p53_mut, size=!is.na(mut), fill=meth_value), color="black") +
        ggrepel::geom_text_repel(aes(label=Name), size=1, alpha=0.5, segment.alpha=0.2, max.overlaps=Inf) +
        facet_grid(meth_id ~ cohort, scales="free") +
        scale_fill_distiller(palette="RdBu", direction=-1) +
#        guides(fill = guide_legend(override.aes=list(shape=21, size=5))) +
        scale_shape_manual(name="p53 mutation", guide="legend", na.value=21,
                           values=mut_shapes) +
        scale_size_manual(name="gene mutation", values=setNames(c(1,3), c("TRUE", "FALSE"))) +
        labs(title = paste("CCLE compensation;",
                           "98/99th% shown (expr/copies); yellow=euploid"),
             y = "normalized read count") +
        theme_classic()
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/meth-ccle.pdf')
    )

    cohorts = c("BRCA", "COADREAD", "LUAD", "LUSC", "SKCM")
    cd = util$load_ccle(c(args$gene,"TP53")) %>%
        filter(gene == args$gene,
               cohort %in% cohorts) %>%
        select(-meth)

    cpg_tss = readr::read_tsv("../data/ccle/CCLE_RRBS_TSS_1kb_20180614.txt") %>%
        filter(gene == args$gene) %>%
        tidyr::gather("CCLE_ID", "meth_value", -(TSS_id:avg_coverage))
    cpg_clust = readr::read_tsv("../data/ccle/CCLE_RRBS_TSS_CpG_clusters_20180614.txt") %>%
        filter(gene_name == args$gene) %>%
        tidyr::gather("CCLE_ID", "meth_value", -(cluster_id:avg_coverage))

    both = bind_rows(cpg_tss %>% select(CCLE_ID, meth_id=TSS_id, meth_value),
                     cpg_clust %>% select(CCLE_ID, meth_id=cluster_id, meth_value))
    comb = inner_join(cd, both)

    pdf(args$plotfile, 16, 15)
    print(plot_ccle_meth(comb))
    dev.off()
})
