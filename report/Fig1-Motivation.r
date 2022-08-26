library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga')
cm = import('./common')

hlg = c("MYC", "EGFR", "CCND1", "CDKN1A", "TP53", "BAP1", "CDKN1A", "IL7R", "CKS1B",
        "APC", "CDKN2A", "KRAS", "NRAS", "RB1", "CCNE1", "PIK3CA", "AURKA")

schema = function() {
    img = grid::rasterGrob(magick::image_read("external/comp+tox.svg"))
    ggplot() + annotation_custom(img) + theme(panel.background=element_blank())
}

overlap = function() {
    img = grid::rasterGrob(magick::image_read("external/overlap.svg"))
    ggplot() + annotation_custom(img) + theme(panel.background=element_blank())
}

cna_along_genome = function(gistic, hlg=c()) {
    labs = gistic$genes %>% filter(gene_name %in% hlg) %>%
        group_by(gene_name) %>%
            slice_max(abs(frac)) %>%
        ungroup() %>%
        inner_join(gistic$smooth %>% select(type, chr, gam)) %>%
        rowwise() %>%
        mutate(frac = mgcv::predict.gam(gam, newdata=data.frame(tss=tss)))

    smooth = gistic$smooth %>% select(-gam) %>% tidyr::unnest(steps) %>%
        mutate(frac_amp = ifelse(frac > 0.15, frac, NA),
               type = stringr::str_to_title(type))
    ggplot(smooth, aes(x=tss)) +
        geom_hline(yintercept=0, color="black") +
        geom_hline(yintercept=0.15, color="black", linetype="dashed") +
        geom_area(aes(y=frac, group=type, fill=type), alpha=0.5) +
        scale_fill_manual(values=cm$cols[c("Amplification", "Deletion")], name="CNA") +
        geom_line(aes(y=frac_amp, group=type, color="Frequently\namplified"),
                  lineend="round", size=1) +
        scale_color_manual(values=c("Frequently\namplified"="#960019"), name="") +
        geom_point(data=labs, aes(y=frac), color="black", fill="white", shape=21) +
        ggrepel::geom_text_repel(data=labs, aes(y=frac, label=gene_name), size=3,
                                 point.size=10, max.iter=1e5, max.time=10) +
        facet_grid(. ~ chr, scales="free", space="free") +
        labs(y = "Alteration frequency") +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.spacing.x = unit(1, "mm")) +
        coord_cartesian(clip="off", expand=FALSE)
}

get_cosmic_annot = function() {
    manual = readRDS("../data/genesets/manual.rds")
    manual[grepl("Cosmic_(OG|TSG)_Tier", names(manual))] %>%
        stack() %>% as_tibble() %>%
        transmute(gene_name = values,
                  type = case_when(
                      grepl("OG", ind) ~ "Oncogene",
                      grepl("TSG", ind) ~ "TSG"
                  ),
                  tier = sub(".*Tier([12])$", "\\1", ind)) %>%
        distinct() %>%
        group_by(gene_name, tier) %>%
        summarize(type = ifelse(length(type) == 1, as.character(type), "Both"))
}

sys$run({
    gistic = readRDS("../data/gistic_smooth.rds")

    top = (schema() | overlap()) + plot_layout(widths=c(3,2))
    btm = wrap_elements(cna_along_genome(gistic, hlg))

    asm = (btm / top) + plot_layout(heights=c(2,3)) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("Fig1-Motivation.pdf", 14, 8)
    print(asm)
    dev.off()
})
