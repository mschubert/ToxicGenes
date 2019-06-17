library(cowplot)
io = import('io')

tissues = io$read_table("tissues.txt", header=TRUE)

plot_overview = function(data, title, label_top=20) {
    data %>%
        mutate(label = ifelse(rank(`LFC DMSO/ETP`) <= label_top, `GENE SYMBOL`, NA)) %>%
        ggplot(aes(x=DMSO, y=`LFC DMSO/ETP`)) +
            geom_point(alpha=0.5) +
            ggrepel::geom_label_repel(aes(label=label)) +
            labs(title = title,
                 x = "DMSO",
                 y = "LFC DMSO/ETP")
}

expr = io$read_table("./data/ORF_DMSO_2019-02.txt", header=TRUE) %>%
    tidyr::gather("condition", "value", -(`Construct Barcode`:`BEST GENE MATCH`)) %>%
    mutate(condition = sub("( ORF)?[_-]DMSO", " DMSO", condition),
           condition = sub("LFC$", "LFC DMSO/ETP", condition),
           cells = sub(" .*$", "", condition),
           cells = sub("LnCAP", "LnCaP", cells, fixed=TRUE),
           cells = sub("SKNEP1", "SK-NEP-1", cells, fixed=TRUE),
           condition = sub("^[A-Za-z0-9-]+ ", "", condition)) %>%
    tidyr::spread(condition, value) %>%
    left_join(tissues %>% select(-comment))

percell = expr %>%
    group_by(cells) %>%
    tidyr::nest() %>%
    mutate(plot = purrr::map2(data, cells, plot_overview))

pdf("overview.pdf")
for (p in percell$plot)
    print(p)
dev.off()
