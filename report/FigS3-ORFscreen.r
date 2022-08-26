library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
orf = import('../orf/overview_naive')

facet_plot = function(ov) {
    ggplot(ov, aes(x=DMSO, y=`LFC DMSO/ETP`)) +
        geom_density2d() +
        orf$stat_loess_sd(color="red", size=1) +
#        facet_wrap(~ cells) +
        theme_minimal()
}

plot_norm = function(percell) {
    plots_naive = purrr::map2(percell$data, percell$cells, orf$plot_overview)
    plots_corr =  purrr::map2(percell$data, percell$cells, orf$plot_overview_corrected)
    wrap_plots(plots_naive, ncol=6, tag_level="new") /
        wrap_plots(plots_corr, ncol=6, tag_level="new")
}

# cor between screens

# volc GO

# genes covered by ORF vs. genome

sys$run({
    ov = readRDS("../orf/overview.rds") %>%
        mutate(cells = sprintf("%s (%s)", cells, tissue))
    percell = ov %>%
        group_by(cells) %>%
        tidyr::nest()

#    orfdata = readxl::read_xlsx("../orf/fits_naive.xlsx", sheet="pan") %>%
#        dplyr::rename(gene_name = `GENE SYMBOL`) %>%
#        filter(gene_name != "LOC254896") # not in tcga/ccle data

    norms = plot_norm(percell)

    asm = norms + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS3-ORFscreen.pdf", 14, 16)
    print(asm)
    dev.off()
})
