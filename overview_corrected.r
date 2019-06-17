library(cowplot)
io = import('io')
ov = import('./overview')

plot_overview = function(data, title, label_top=20) {
    data = data %>%
        mutate(label = ifelse(rank(z_LFC) <= label_top, `GENE SYMBOL`, NA))
    ggplot(data, aes(x=DMSO, y=z_LFC)) +
        geom_point(alpha=0.5) +
        ov$loess_sd(data$DMSO, data$z_LFC) +
#        stat_loess_sd(color="red", size=1) +
        ggrepel::geom_label_repel(aes(label=label), size=2) +
        labs(title = title,
             x = "DMSO",
             y = "LFC DMSO/ETP z-score")
}

expr = readRDS("overview.rds")

percell = expr %>%
    group_by(cells) %>%
    tidyr::nest() %>%
    mutate(plot = purrr::map2(data, cells, plot_overview))

pdf("overview_corrected.pdf")
for (p in percell$plot)
    print(p)
dev.off()
