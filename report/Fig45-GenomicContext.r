library(dplyr)
library(plyranges)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
cm = import('./common')

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

    # -> add ARGOS to type (+incl in label repel)
    # -> add events spanning CDKN1A
    labs = gistic$genes %>%
        inner_join(cosmic) %>%
        filter(chr == "6",
               (gtype == "Oncogene" & type == "amplification") |
                (gtype == "TSG" & type == "deletion") |
                gtype == "OG+TSG") %>%
        inner_join(gistic$smooth %>% select(type, chr, gam)) %>%
        rowwise() %>%
        mutate(frac = mgcv::predict.gam(gam, newdata=data.frame(tss=tss)))
    labs2 = labs %>% filter(gene_name != "CDKN1A")
    labs3 = labs %>% filter(gene_name == "CDKN1A") %>%
        select(gene_name, type, frac, tss) %>%
        tidyr::spread(type, frac)

    p21 = gistic$smooth %>%
        filter(chr == "6") %>%
        select(-gam) %>% tidyr::unnest(steps) %>%
        mutate(frac_amp = ifelse(frac > 0.15, frac, NA),
               type = stringr::str_to_title(type))

    ev21 = join_overlap_inner(ev, genes %>% filter(gene_name == "CDKN1A")) %>%
        as.data.frame() %>% as_tibble() %>%
        filter(event_type == "amp") %>%
        mutate(frac_start = rank(start, ties.method="first")/nrow(.),
               frac_end = rank(-end, ties.method="last")/nrow(.))
    pev = ggplot(ev21) +
        geom_vline(xintercept=labs3$tss, color=cm$cols["Compensated"], linewidth=2, alpha=0.5) +
        geom_step(aes(x=start, y=frac_start)) +
        geom_step(aes(x=end, y=frac_end)) +
        annotate("point", x=labs3$tss, y=1, size=5, shape=21,
                 fill=cm$cols["Compensated"], alpha=0.9) +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank()) +
        ylab("Event length") +
        coord_cartesian(clip="off")

    pcn = ggplot(p21, aes(x=tss)) +
        geom_segment(data=labs3, aes(y=amplification, yend=deletion, x=tss, xend=tss),
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
        scale_fill_manual(values=cm$cols, name="Driver") +
        ggrepel::geom_text_repel(data=labs2, aes(y=frac, label=gene_name), size=3,
            point.size=5, max.iter=1e5, max.time=10, segment.alpha=0.3) +
        annotate("point", x=labs3$tss, y=0, size=5, shape=21,
                 fill=cm$cols["Compensated"], alpha=0.9) +
        facet_grid(. ~ chr, scales="free", space="free") +
        labs(y = "Alteration frequency TCGA") +
        theme_minimal() +
        theme(axis.title.x = element_blank())

    (pcn/pev) + plot_layout(heights=c(3,1), guides="collect")
})
