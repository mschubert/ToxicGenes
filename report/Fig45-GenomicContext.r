library(dplyr)
library(plyranges)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
cm = import('./common')

plot_ctx = function(genes, ev, cosmic, gistic, .hl) {
    cur_ev = join_overlap_inner(ev, genes %>% filter(gene_name == .hl)) %>%
        as.data.frame() %>% as_tibble() %>%
        filter(event_type == "amp") %>%
        mutate(frac_start = rank(start, ties.method="first")/nrow(.),
               frac_end = rank(-end, ties.method="last")/nrow(.))

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
    labs2 = labs %>% filter(!duplicated(gene_name)) %>%
        mutate(frac = ifelse(gene_name == .hl, 0, frac),
               gtype = ifelse(gene_name == .hl, "ARGOS", gtype))
    rng = labs %>% filter(gene_name == .hl) %>%
        select(gene_name, type, frac, tss) %>%
        tidyr::spread(type, frac)

    cnv = gistic$smooth %>%
        filter(chr == cur_ev$seqnames[1]) %>%
        select(-gam) %>% tidyr::unnest(steps) %>%
        mutate(frac_amp = ifelse(frac > 0.15, frac, NA),
               type = stringr::str_to_title(type))

    pev = ggplot(cur_ev) +
        geom_vline(xintercept=rng$tss, color=cm$cols["Compensated"], linewidth=2, alpha=0.5) +
        geom_step(aes(x=start, y=frac_start)) +
        geom_step(aes(x=end, y=frac_end)) +
        annotate("point", x=rng$tss, y=1, size=5, shape=21,
                 fill=cm$cols["Compensated"], alpha=0.9) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank()) +
        ylab("Event fraction") +
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
        geom_point(data=labs2, aes(y=frac, fill=gtype), alpha=0.8, color="black", shape=21, size=5) +
        scale_fill_manual(values=c(cm$cols, ARGOS=cm$cols[["Compensated"]]), name="Driver") +
        ggrepel::geom_text_repel(data=labs2, aes(y=frac, label=gene_name), size=3,
            point.size=5, max.iter=1e5, max.time=10, segment.alpha=0.3) +
        annotate("point", x=rng$tss, y=0, size=5, shape=21,
                 fill=cm$cols["Compensated"], alpha=0.9) +
        facet_grid(. ~ chr, scales="free", space="free") +
        labs(y = "Alteration frequency TCGA") +
        theme_minimal() +
        theme(axis.title.x = element_blank())

    (pcn/pev) + plot_layout(heights=c(3,1), guides="collect")

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
        filter(abs(cn_end - cn_start) > 0.5) %>%
        transmute(Sample=sample, seqnames=chr, start=base_start, end=base_end,
                  seg_id = seq_len(nrow(.)), event_type=event_type) %>%
        makeGRangesFromDataFrame(keep.extra.columns=TRUE)
    gistic = readRDS("../data/gistic_smooth.rds")

    pdf("Fig45-GenomicContext.pdf", 6, 4)
    print(plot_ctx(genes, ev, cosmic, gistic, "CDKN1A"))
    print(plot_ctx(genes, ev, cosmic, gistic, "RBM14"))
    dev.off()
})
