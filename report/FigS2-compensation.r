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

volcs = function() {
    ccle_genes = readRDS("../ccle/pan/stan-nb.rds")$amp
    tcga_genes = readRDS("../tcga/pan/stan-nb_puradj.rds")
    ccle_go = readRDS("../ccle/pan/stan-nb/all/GO_Biological_Process_2021.rds")$amp
    tcga_go = readRDS("../tcga/pan/stan-nb_puradj/all/GO_Biological_Process_2021.rds")[[1]]

    padj = "Adjusted p-values (FDR)"
    ccle = (plt$volcano(ccle_genes) + guides(size="none")) |
        (plt$volcano(ccle_go) + labs(size="Set size", y=padj)) +
        plot_layout(guides="collect")
    tcga = (plt$volcano(tcga_genes) + guides(size="none")) |
        (plt$volcano(tcga_go, pos_label_bias=0.5) + labs(size="Set size", y=padj)) +
        plot_layout(guides="collect")

    ccle / tcga
}

# mcmc traces of some example genes

# + cor plot tcga vs ccle for gene sets

sys$run({
    asm = (tcga_vs_ccle() / volcs()) + plot_layout(heights=c(1,4))

    asm = asm + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("FigS2-compensation.pdf", 10, 14)
    print(asm)
    dev.off()
})
