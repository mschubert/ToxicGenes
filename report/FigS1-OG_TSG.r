library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
cm = import('./common')

og_tsg_cna = function(gistic, cosmic) {
    dset = gistic %>%
        mutate(frac = abs(frac),
               cna = stringr::str_to_title(type)) %>% select(-type) %>%
        left_join(cosmic) %>%
        mutate(type = ifelse(is.na(type), "Background", type),
               type = factor(type, levels=c("Background", "Oncogene", "TSG"))) %>%
        filter(!is.na(type))

    bg_line = dset %>% group_by(cna) %>%
        summarize(frac = median(frac[type == "Background"], na.rm=TRUE))

    ggplot(dset, aes(x=type, y=frac, fill=type)) +
        geom_boxplot(outlier.shape=NA) +
        scale_fill_manual(values=cm$cols[levels(dset$type)]) +
        labs(y="Frequency TCGA", x ="Driver status (whole genome)",
             fill="Driver status\n(whole genome)") +
        geom_hline(data=bg_line, aes(yintercept=frac), linetype="dashed", color="black") +
        ggsignif::geom_signif(comparisons=list(c("Background", "Oncogene"), c("Background", "TSG")),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            y_position=c(0.43,0.48,0.43,0.48), color="black", test=t.test) +
        facet_wrap(~ cna) +
        coord_cartesian(ylim=c(0.02, 0.52)) +
        cm$theme_classic() +
        theme(strip.background = element_blank(),
              strip.placement = "outside",
              axis.text.x = element_blank())
}

og_vs_tsg = function(gistic, cosmic) {
    gwide = tidyr::pivot_wider(gistic, names_from="type", values_from="frac") %>%
        left_join(cosmic) %>%
        filter(!is.na(type)) %>%
        mutate(label = ifelse(!is.na(type) & gene_name %in% cm$hlg, gene_name, NA))

    areas = data.frame(xmin=c(0.15,-Inf), xmax=c(Inf,Inf), ymin=c(-Inf,0.15), ymax=c(Inf,Inf),
        status=c("Deleted", "Amplified"))

    ggplot(gwide, aes(x=-deletion, y=amplification, color=type)) +
        geom_rect(data=areas, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=status),
            linetype="dashed", alpha=0.1, inherit.aes=FALSE) +
        scale_fill_manual(values=cm$cols[c("Amplified", "Deleted")], name="Frequently") +
        geom_point(aes(shape=tier), na.rm=TRUE) +
        ggrepel::geom_label_repel(aes(label=label), min.segment.length=0,
            segment.alpha=0.3, fill="#ffffffca", label.size=NA, na.rm=TRUE) +
        scale_shape_manual(values=c("1"=19, "2"=1), name="COSMIC tier") +
        scale_color_manual(values=cm$cols[c("Oncogene", "TSG", "OG+TSG")], name="Driver status") +
        coord_fixed() +
        cm$theme_classic() +
        labs(x = "Deletion frequency",
             y = "Amplification frequency")
}

venn_amp_del = function(gistic, cosmic) {
    amps = gistic %>% filter(type == "amplification")
    dels =  gistic %>% filter(type == "deletion")
    sets = list(
        `Frequently amplified` = unique(amps$gene_name[amps$frac > 0.15]),
        `Frequently deleted` = unique(dels$gene_name[dels$frac < -0.15])
    )
    sets$Neither = setdiff(unique(gistic$gene_name), unlist(sets, use.names=FALSE))
    cols = cm$cols[c("Amplified", "Deleted", "Background")]
    names(cols) = names(sets)
    plt$venn(sets, alpha=0.3) +
        scale_fill_manual(values=cols) +
        theme(plot.margin = margin(0,5,0,5, unit="mm"))
}

venn_og_tsg = function(gistic, cosmic) {
    amps = gistic %>% filter(type == "amplification")
    dels =  gistic %>% filter(type == "deletion")
    sets = list(
        `Frequently\namplified` = unique(amps$gene_name[amps$frac > 0.15]),
        Oncogene = cosmic$gene_name[grepl("^O", cosmic$type)],
        TSG = cosmic$gene_name[grepl("TSG", cosmic$type)]
    )
    cols = cm$cols[c("Oncogene", "TSG", "Amplified")]
    names(cols)[3] = names(sets)[1]
    plt$venn(sets, alpha=0.3) +
        scale_fill_manual(values=cols) +
        theme(plot.margin = margin(0,5,0,5, unit="mm"))
}

rpe_scaling = function(rpe) {
    lookup = c("14.10"="14.10 (+7 +16 +X)", "14.16"="14.16 (+20)", "14.21"="14.21 (+8)")
    diff_expr = rpe$diff_expr %>% mutate(clone = lookup[clone])
    segs = rpe$segs %>% mutate(clone = lookup[clone])

    p1 = ggplot(diff_expr, aes(x=loc, y=log2FoldChange)) +
        geom_hline(yintercept=c(log2((1:4)/2)), color="firebrick", linetype="dashed") +
        geom_point(size=0.5, alpha=0.1) +
        geom_segment(data=segs, aes(color=type, y=lfc, yend=lfc),
                     x=-Inf, xend=Inf, size=1.5, alpha=0.8) +
        facet_grid(clone ~ seqnames, scales="free", space="free") +
        coord_cartesian(ylim=c(-1.8,1.8)) +
        scale_color_manual(values=c(DNA="purple", RNA="green"), name="Data type") +
        cm$theme_minimal() +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.background = element_rect(fill="#f5f5f5"),
              panel.spacing.x = unit(1, "mm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +
        xlab("Genomic location")

    comp = segs %>% group_by(clone, seqnames, seg_id) %>%
        summarize(is_amp = cut(lfc[type=="DNA"], c(-Inf, -0.25, 0.25, Inf),
                               labels=c("Deleted", "Euploid", "Amplified")),
                  lfc_diff = lfc[type=="RNA"]-lfc[type=="DNA"]) %>%
        filter(is_amp != "Deleted", seqnames != "10")

    p2 = ggplot(comp, aes(x=is_amp, y=lfc_diff)) +
        geom_boxplot(aes(fill=is_amp), outlier.shape=NA, alpha=0.5) +
        scale_fill_manual(values=cm$cols[c("Euploid", "Amplified")], guide="none") +
        scale_color_brewer(palette="Dark2") +
        ggbeeswarm::geom_quasirandom(size=2, aes(color=clone)) +
        ggsignif::geom_signif(comparisons=list(c("Euploid", "Amplified")),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            y_position=c(0.25,0.28), color="black", test=t.test) +
        cm$theme_classic() +
        coord_cartesian(clip="off") +
        ylab("LFC RNA/DNA chromosome") +
        theme(axis.title.x = element_blank())

    list(genome=p1, quant=p2)
}

# cleanup: copied/adapted from dorine/data/dna_vs_rna.r
center_maxdens = function(y, bw=bw.nrd0) {
    max_dens = function(y, bw=bw.nrd0) {
        dens = density(y, kernel="gaussian")
        dens$x[dens$y==max(dens$y)]
    }
    y - max_dens(y, bw=bw)
}

rpe2_scaling = function() {
    lookup = c(SS6="SS6 (+7)", SS51="SS51 (+7 +22)", SS111="SS111 (+8 +9 +18)")
    rna = readRDS("../data/rnaseq_rpe1_broad/compute_fcs.rds") %>%
        tidyr::unnest(genes) %>%
        filter(lfcSE < 4) %>%
        transmute(Sample=lookup[term], Gene=label, FC=log2FoldChange, lfcSE=lfcSE)
#    means = function(mat) narray::map(mat, along=2, mean, subsets=sub("-[0-9]+$", "", colnames(mat)))
#    rpe2 = readxl::read_xlsx("../data/Expression-matrix_RPE1-clones_reads.xlsx", skip=1)
#    rna = rpe2 %>%
#        transmute(Gene=Gene, chr=factor(sub("[pq].*$", "", Location), levels=c(1:22,'X'))) %>%
#        cbind(means(data.matrix(rpe2[-c(1,2)]))) %>% as_tibble() %>%
#        tidyr::gather("Sample", "expr", -Gene, -chr, -SS48) %>%
#        mutate(Sample = factor(lookup[Sample], levels=c(lookup)),
#               FC = expr / SS48) %>%
#        filter(SS48 >= 20 & expr >= 20, !is.na(Sample))
    dna = readxl::read_xlsx("../data/Supp Table 3 - Gene Copy Number analysis.xlsx") %>%
        mutate(chr=factor(sub("[pq].*$", "", `Chromosome localization`), levels=c(1:22,'X'))) %>%
        select(Gene = `Gene name`, chr, loc=`Gene start (bp)`, `RPE1-SS48`:`RPE1-SS111`) %>%
        tidyr::gather("Sample", "copy", -Gene, -chr, -loc) %>%
        mutate(Sample = lookup[sub("^RPE1-", "", Sample)]) %>%
        na.omit()
    both = inner_join(rna, dna) %>% filter(chr != "10")
    dna_segs = both %>% group_by(Sample, chr) %>% summarize(mean_copy=mean(copy))
    rna_segs = both %>% group_by(Sample, chr) %>%
        summarize(mean_FC=mean(FC, weights=1/lfcSE, trim=0.2)) %>%
        group_by(Sample) %>% mutate(mean_FC = center_maxdens(mean_FC))

    rna_plot = ggplot(both, aes(x=loc, y=2^FC)) + geom_point(alpha=0.1) +
        facet_grid(Sample~chr, scale="free", space="free") + scale_y_log10() +
        geom_hline(yintercept=c(1,2,3,4)/2, color="firebrick", linetype="dashed") +
        geom_hline(data=dna_segs, aes(yintercept=mean_copy), color="purple", linewidth=1) +
        geom_hline(data=rna_segs, aes(yintercept=2^mean_FC), color="green", linewidth=1) +
        coord_cartesian(ylim=c(0.15,5)) + ggtitle("RNA")
    dna_plot = ggplot(both, aes(x=loc, y=copy)) + geom_point(alpha=0.1) +
        facet_grid(Sample~chr, scale="free", space="free") + scale_y_log10() +
        geom_hline(yintercept=c(1,2,3,4)/2, color="firebrick", linetype="dashed") +
        geom_hline(data=dna_segs, aes(yintercept=mean_copy), color="purple", linewidth=1) +
        coord_cartesian(ylim=c(0.15,5)) + ggtitle("DNA")

    comp = inner_join(rna_segs, dna_segs) %>%
        mutate(is_amp = ifelse(mean_copy<1.2, "Euploid", "Amplified") %>%
                    factor(levels=c("Euploid", "Amplified")),
               lfc_diff = mean_FC-log2(mean_copy))
    quant_plot = ggplot(comp, aes(x=is_amp, y=lfc_diff)) +
        geom_boxplot(aes(fill=is_amp), outlier.shape=NA, alpha=0.5) +
        scale_fill_manual(values=cm$cols[c("Euploid", "Amplified")], guide="none") +
        scale_color_brewer(palette="Dark2") +
        ggbeeswarm::geom_quasirandom(size=2, aes(color=Sample)) +
        ggsignif::geom_signif(comparisons=list(c("Euploid", "Amplified")),
            map_signif_level=cm$fmt_p, parse=TRUE, tip_length=0,
            y_position=c(0.25,0.28), color="black", test=t.test, textsize=3) +
        cm$theme_classic() +
        coord_cartesian(clip="off") +
        ylab("LFC RNA/DNA chromosome") +
        theme(axis.title.x = element_blank())

    list(genome=rna_plot, quant=quant_plot)
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")$genes
    cosmic = cm$get_cosmic_annot()
    rpe = readRDS("../data/dorine_compare.rds")
    rs = rpe_scaling(rpe)

    left = (og_vs_tsg(gistic, cosmic) / og_tsg_cna(gistic, cosmic)) +
        plot_layout(heights=c(2,1))

#    right = (venn_amp_del(gistic, cosmic) / venn_og_tsg(gistic, cosmic) / rs$quant) +
    right = (venn_amp_del(gistic, cosmic) / venn_og_tsg(gistic, cosmic)) +
        plot_layout(heights=c(3,4))

#    asm = ((((left | right) + plot_layout(widths=c(3,2))) / rs$genome) + plot_layout(heights=c(3,2))) +
    asm = ((left | right) + plot_layout(widths=c(3,2))) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=24, face="bold"))

    pdf("FigS1-OG_TSG.pdf", 11, 8)
    print(asm)
    dev.off()
})
