library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

tcga_vs_ccle = function() {
    ccle = readxl::read_xlsx("../ccle/pan/stan-nb.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga1 = readxl::read_xlsx("../tcga/pan/stan-nb_naive.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga2 = readxl::read_xlsx("../tcga/pan/stan-nb_pur.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))
    tcga3 = readxl::read_xlsx("../tcga/pan/stan-nb_puradj.xlsx") %>%
        mutate(estimate = pmax(-2, pmin((1 - p.value) * estimate, 2.5)))

    dset = ccle %>% select(gene, CCLE=estimate) %>%
        left_join(tcga1 %>% select(gene, `No purity correction`=estimate)) %>%
        left_join(tcga2 %>% select(gene, `Common purity term`=estimate)) %>%
        left_join(tcga3 %>% select(gene, `Purity per tissue`=estimate)) %>%
        tidyr::gather("type", "value", -gene, -CCLE) %>%
        mutate(type = factor(type) %>% relevel("No purity correction"))

    mods = dset %>% group_by(type) %>%
        summarize(mod = list(lm(value ~ CCLE))) %>%
        rowwise() %>%
        mutate(tidy = list(broom::tidy(mod)),
               glance = list(broom::glance(mod)),
               intcp = tidy$estimate[tidy$term == "(Intercept)"],
               slope = tidy$estimate[tidy$term == "CCLE"],
               angle = atan(slope) * 180/pi) %>%
        select(-tidy, -mod) %>%
        tidyr::unnest(glance) %>%
        mutate(label = sprintf("R^2~`=`~%.2f~p~`=`~10^%.0f", adj.r.squared, ceiling(log10(p.value))))

    ggplot(dset, aes(x=CCLE, y=value)) +
        geom_vline(xintercept=0, color="grey", linetype="dashed", size=1) +
        geom_hline(yintercept=0, color="grey", linetype="dashed", size=1) +
        geom_hex(aes(color=..count..), bins=50) +
        scale_color_continuous(type = "viridis", trans="log1p", guide="none") +
        scale_fill_continuous(type = "viridis", trans="log1p", breaks=c(1,5,20,100,500)) +
        facet_wrap(~ type) +
        geom_smooth(method="lm", color="red", se=FALSE, size=0.7) +
        geom_text(data=mods, aes(x=0, y=intcp, label=label, angle=angle), parse=TRUE,
                  color="red", hjust=0.4, vjust=-0.5, size=3) +
        labs(y = "Expression over expected TCGA",
             x = "Expression over expected CCLE") +
        coord_cartesian(ylim=c(-1.2, 1.2)) +
        theme_minimal() +
        theme(strip.text = element_text(size=12),
              strip.background = element_rect(color=NA, fill="#ffffffc0"))
}

go_cors = function() {
    ccle_go = readRDS("../ccle/pan/stan-nb/all/GO_Biological_Process_2021.rds")$amp %>%
        select(label, stat_ccle=statistic)
    tcga_go = readRDS("../tcga/pan/stan-nb_puradj/all/GO_Biological_Process_2021.rds")[[1]] %>%
        select(label, stat_tcga=statistic, size_used)
    both = inner_join(ccle_go, tcga_go) %>% filter(size_used < 1000)

    m = broom::glance(lm(stat_ccle ~ stat_tcga, data=both))
    lab = sprintf("R^2~`=`~%.2f~\n~p~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    plt$denspt(both, aes(x=stat_tcga, y=stat_ccle, label=label), size=size_used) +
        scale_size_area(max_size=8, breaks=c(10,100,500,1000), name="Genes in set") +
        theme_minimal() +
        labs(title = "Gene Ontology: Biological Process",
             x = "Δ Expression over expected TCGA (Wald stat.)",
             y = "Δ Expression over expected CCLE (Wald stat.)") +
        annotate("text", x=11, y=-6.5, color="blue", label=lab, parse=TRUE)
}

rpe_comp = function(rpe, all) {
    gclass = all %>%
        dplyr::rename(label = gene) %>%
        mutate(gclass = case_when(
            est_ccle < -0.3 & est_tcga < -0.3 ~ "Compensated",
            est_ccle > 0.3 & est_tcga > 0.3 ~ "Hyperactivated",
#            abs(est_ccle) < 0.3 & abs(est_tcga) < 0.3 ~ "Background"
            TRUE ~ "Background"
        ))

    comp2 = rpe$segs %>% filter(type == "DNA") %>%
        inner_join(rpe$diff_expr, by=c("clone", "seqnames")) %>%
        mutate(cna = cut(lfc[type=="DNA"], c(-Inf, -0.15, 0.15, Inf),
                            labels=c("Deleted", "Euploid", "Amplified")),
               lfc_diff = log2FoldChange-lfc) %>%
        group_by(seqnames) %>%
            mutate(chr_has_amp = any(cna == "Amplified")) %>%
        ungroup() %>%
        inner_join(gclass) %>%
        mutate(group = case_when(
            chr_has_amp & cna == "Euploid" & gclass == "Background" ~ "Euploid\nchr 8,12,13,16,20",
            cna == "Euploid" & gclass == "Background" ~ "Background\nother chr",
            cna == "Amplified" & gclass == "Background" ~ "Amplified\nNon-Comp.",
            cna == "Amplified" & gclass == "Compensated" ~ "Amplified\nCompensated"
        )) %>% filter(!is.na(group)) %>%
            mutate(group = factor(group, levels=c("Background\nother chr",
                "Euploid\nchr 8,12,13,16,20", "Amplified\nNon-Comp.", "Amplified\nCompensated")))

    ggplot(comp2, aes(x=group, y=lfc_diff)) +
        geom_boxplot(outlier.shape=NA) +
        geom_hline(yintercept=0, linetype="dashed", color="grey") +
        coord_cartesian(ylim=c(-2,2.8), clip="off") +
        theme_classic() +
        ggsignif::geom_signif(comparisons=list(
                c("Background\nother chr", "Euploid\nchr 8,12,13,16,20"),
                c("Background\nother chr", "Amplified\nNon-Comp."),
                c("Background\nother chr", "Amplified\nCompensated"),
                c("Amplified\nNon-Comp.", "Amplified\nCompensated")),
            y_position=c(1.8,1.5,1.2,0.9), color="black", test=t.test, textsize=3,
            tip_length=0.002) +
        labs(x = "Group",
             y = "LFC DNA/RNA isogenic RPE-1 clones")
}

# mcmc traces of some example genes

# + cor plot tcga vs ccle for gene sets

sys$run({
    rpe = readRDS("../data/dorine_compare.rds")
    all = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")

    left = (tcga_vs_ccle() / go_cors()) + plot_layout(heights=c(1,3))
    right = (plot_spacer() / rpe_comp(rpe, all)) + plot_layout(heights=c(2,1))

    asm = (left | right) + plot_layout(widths=c(2,1)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("FigS2-Compensation.pdf", 15, 12)
    print(asm)
    dev.off()
})
