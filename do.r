library(cowplot)
io = import('io')

plot_overview = function(data, label_top=20) {
    data %>%
        mutate(label = ifelse(rank(`LFC DMSO/ETP`) <= label_top, `GENE SYMBOL`, NA)) %>%
        ggplot(aes(x=DMSO, y=`LFC DMSO/ETP`)) +
            geom_point(alpha=0.5) +
            ggrepel::geom_label_repel(aes(label=label))
}

expr = io$read_table("./data/ORF_DMSO_2019-02.txt", header=TRUE) %>%
    tidyr::gather("condition", "value", -(`Construct Barcode`:`BEST GENE MATCH`)) %>%
    mutate(condition = sub("( ORF)?[_-]DMSO", " DMSO", condition),
           condition = sub("LFC$", "LFC DMSO/ETP", condition),
           cells = sub(" .*$", "", condition),
           condition = sub("^[A-Za-z0-9-]+ ", "", condition)) %>%
    tidyr::spread(condition, value)

percell = expr %>%
    group_by(cells) %>%
    tidyr::nest() %>%
    purrr::map(data, plot_overview)
