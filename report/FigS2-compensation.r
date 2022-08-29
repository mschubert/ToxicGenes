library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')

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
        left_join(tcga2 %>% select(gene, `Purity overall`=estimate)) %>%
        left_join(tcga3 %>% select(gene, `Purity per tissue`=estimate)) %>%
        tidyr::gather("type", "value", -gene, -CCLE)

    ggplot(dset, aes(x=CCLE, y=value)) +
        geom_hex(aes(color=..count..), bins=50) +
        scale_color_continuous(type = "viridis", trans="log1p", guide="none") +
        scale_fill_continuous(type = "viridis", trans="log1p", breaks=c(1,5,20,100,500)) +
        facet_wrap(~ type) +
        geom_smooth(method="lm", color="red", se=FALSE) +
        labs(x = "Expression over expected TCGA",
             y = "Expression over expected CCLE") +
        coord_cartesian(ylim=c(-1.2, 1.2)) +
        theme_minimal() +
        theme(strip.background = element_rect(color=NA, fill="#ffffffc0"))
}

# tcga-vs-ccle comp purity correction plots

# mcmc traces of some example genes

# volcano plots for GSEA (+genes?)
# + cor plot tcga vs ccle for gene sets

sys$run({
    asm = tcga_vs_ccle()

    asm = asm + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    cairo_pdf("FigS2-compensation.pdf", 14, 5)
    print(asm)
    dev.off()
})
