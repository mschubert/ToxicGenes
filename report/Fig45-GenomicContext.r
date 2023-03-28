library(dplyr)
library(plyranges)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
cm = import('./common')
tcga = import('data/tcga')

tissue_compare = function(.hl) {
    # todo: fix column names in result files, save rds as well
    load_orf = . %>% readxl::read_xlsx() %>% dplyr::rename(gene=name)
    dset = list(
        CCLE = list(`Pan-cancer` = readRDS("../ccle/pan/stan-nb.rds")$amp,
                    Lung = readRDS("../ccle/NSCLC/stan-nb.rds"),
                    Breast = readRDS("../ccle/BRCA/stan-nb.rds")),
        TCGA = list(`Pan-cancer` = readRDS("../tcga/NSCLC/stan-nb_puradj.rds"),
                    Lung = readRDS("../tcga/NSCLC/stan-nb_puradj.rds"),
                    Breast = readRDS("../tcga/BRCA/stan-nb_puradj.rds")),
        ORF = list(`Pan-cancer` = readxl::read_xlsx("../orf/fits_naive.xlsx") %>%
                        dplyr::rename(gene=`GENE SYMBOL`),
                   Lung = load_orf("../orf/BRCA/genes.xlsx"),
                   Breast = load_orf("../orf/LUAD/genes.xlsx"))
    ) %>% lapply(bind_rows, .id="tissue") %>% bind_rows(.id="dset") %>%
        mutate(std.error = ifelse(is.na(std.error), estimate/z_comp, std.error),
               Tissue = factor(tissue, levels=c("Pan-cancer", "Breast", "Lung")),
               dset = factor(dset, levels=c("TCGA", "CCLE", "ORF"))) %>%
        filter(gene == .hl) %>% select(Tissue, dset, gene, estimate, std.error)

    ggplot(dset, aes(x=estimate, y=Tissue, fill=Tissue)) +
        geom_col() +
        geom_errorbarh(aes(xmin=estimate-std.error, xmax=estimate+std.error),
                       height=0.2, alpha=0.3) +
        geom_vline(xintercept=0) +
        scale_y_discrete(limits=rev) +
        scale_x_continuous(breaks=c(-0.5, -2)) +
        facet_wrap(~ dset, scales="free_x") +
        scale_fill_brewer(palette="Accent", guide="none") +
        xlab("Compensation (score) / ORF dropout (log2 FC)") +
        theme_minimal() +
        theme(axis.title.y = element_blank())
}

plot_ctx = function(genes, ev, cosmic, gistic, .hl) {
    cur_ev = join_overlap_inner(ev, genes %>% filter(gene_name == .hl)) %>%
        as.data.frame() %>% as_tibble() %>%
        filter(event_type == "amp") %>%
        mutate(cohort = tcga$barcode2study(Sample))
    fracs = . %>% mutate(frac_start = rank(start, ties.method="first"),
                         frac_end = rank(-end, ties.method="last"))
    cur_ev = list(`Pan-cancer` = fracs(cur_ev),
                  Breast = fracs(cur_ev %>% filter(cohort == "BRCA")),
                  Lung = fracs(cur_ev %>% filter(cohort %in% c("BRCA", "LUAD", "LUSC")))) %>%
        bind_rows(.id="Tissue") %>%
        mutate(Tissue = factor(Tissue, levels=c("Pan-cancer", "Breast", "Lung")))

    labs = gistic$genes %>%
        left_join(cosmic) %>%
        filter(chr == cur_ev$seqnames[1],
               (gtype == "Oncogene" & type == "amplification") |
                (gtype == "TSG" & type == "deletion") |
                gene_name == .hl) %>%
        inner_join(gistic$smooth %>% select(type, chr, gam)) %>%
        rowwise() %>%
            mutate(frac = mgcv::predict.gam(gam, newdata=data.frame(tss=tss))) %>%
        ungroup()
    labs2 = labs %>%
        mutate(frac = ifelse(gene_name == .hl, 0, frac),
               gtype = ifelse(gene_name == .hl, "ARGOS", gtype)) %>%
        filter(! (gtype == "ARGOS" & type == "deletion"))
    labs3 = labs2 %>%
        mutate(gene_name = ifelse(hallmark == "Yes" | gene_name %in% c(.hl, "CCND3"),
                                  gene_name, ""))
    rng = labs %>% filter(gene_name == .hl) %>%
        select(gene_name, type, frac, tss) %>%
        tidyr::spread(type, frac)

    cnv = gistic$smooth %>%
        filter(chr == cur_ev$seqnames[1]) %>%
        select(-gam) %>% tidyr::unnest(steps) %>%
        mutate(frac_amp = ifelse(frac > 0.15, frac, NA),
               type = stringr::str_to_title(type))

    pev = ggplot(cur_ev, aes(color=Tissue)) +
        geom_vline(xintercept=rng$tss, color="black", alpha=0.5) +
        geom_step(aes(x=start, y=frac_start)) +
        geom_step(aes(x=end, y=frac_end)) +
        theme_minimal() +
        theme(axis.text.x = element_blank()) +
        labs(x = "Genomic location (bp)",
             y = "Amplification\nevents (cum.)") +
        scale_color_brewer(palette="Accent", guide="none") +
        coord_cartesian(clip="off")

    pcn = ggplot(cnv, aes(x=tss)) +
        geom_segment(data=rng, aes(y=amplification, yend=deletion, x=tss, xend=tss),
                     color=cm$cols["Compensated"], linewidth=2) +
        geom_hline(yintercept=0, color="black") +
        geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
        geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
        scale_fill_manual(values=cm$cols[c("Amplification", "Deletion")], name="CNA") +
        geom_line(aes(y=frac_amp, group=type, color="Frequently\namplified"),
                  lineend="round", size=1) +
        scale_color_manual(values=c("Frequently\namplified"="#960019"), name="") +
        ggnewscale::new_scale("fill") +
        geom_point(data=labs2, aes(y=frac, fill=gtype, size=hallmark),
                   alpha=0.8, color="black", shape=21) +
        scale_size_manual(values=c(Yes=3, No=2), name="COSMIC\nHallmark") +
        scale_fill_manual(values=c(cm$cols, ARGOS=cm$cols[["Compensated"]]), name="Driver") +
        guides(fill = guide_legend(override.aes=list(size=2.5))) +
        ggrepel::geom_text_repel(data=labs3, aes(y=frac, label=gene_name),
            size=3, point.size=3, min.segment.length=0, max.iter=1e5,
            max.time=10, segment.alpha=0.3, force_pull=0.01) +
        annotate("point", x=rng$tss, y=0, size=5, shape=21,
                 fill=cm$cols["Compensated"], alpha=0.9) +
        facet_grid(. ~ chr, scales="free", space="free") +
        labs(x = "Genomic location (bp)",
             y = "Alteration frequency TCGA") +
        theme_minimal() #+
#        theme(axis.title.x = element_blank())

#    (pcn/pev/tissue_compare(.hl)) + plot_layout(heights=c(3,1.2,0.8), guides="collect")
    (pcn/tissue_compare(.hl)) + plot_layout(heights=c(3,0.8), guides="collect")
}

sys$run({
    #todo: save genes, ev in ../data/df_tcga_copysegments.r
    genes = seq$gene_table() %>%
        filter(gene_biotype == "protein_coding") %>%
        group_by(ensembl_gene_id) %>%
        summarize(seqnames=unique(chromosome_name), start=min(start_position),
                  stop=max(end_position), strand=c("-","+")[unique(strand/2+1.5)],
                  gene_name=unique(external_gene_name)) %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    cosmic = cm$get_cosmic_annot() %>% dplyr::rename(gtype=type)
    ev = readr::read_tsv("../data/gistic/TCGA.all_cancers.150601.zigg_events.160923.txt") %>%
        filter(abs(cn_end - cn_start) > 0.2) %>%
        transmute(Sample=sample, seqnames=chr, start=base_start, end=base_end,
                  seg_id = seq_len(nrow(.)), event_type=event_type) %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    gistic = readRDS("../data/gistic_smooth.rds")

    pdf("Fig45-GenomicContext.pdf", 6.5, 4.5)
    print(plot_ctx(genes, ev, cosmic, gistic, "CDKN1A"))
    print(plot_ctx(genes, ev, cosmic, gistic, "RBM14"))
    dev.off()
})
