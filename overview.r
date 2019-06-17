library(cowplot)
io = import('io')

tissues = io$read_table("tissues.txt", header=TRUE)

loess_sd = function(fml, data) {
    x = data[[as.character(fml[[3]])]]
    y = data[[as.character(fml[[2]])]]
    mod = msir::loess.sd(x, y, data=data)
    df = data.frame(x = mod$x, y=mod$y, sd=mod$sd)
    df = df[seq(1, nrow(df), length.out=100),]
    list(geom_line(data=df, aes(x=x, y=y), color="red", size=1),
         geom_line(data=df, aes(x=x, y=y+sd), color="red", linetype="dotted", size=1),
         geom_line(data=df, aes(x=x, y=y-sd), color="red", linetype="dotted", size=1))
}

plot_overview = function(data, title, label_top=20) {
    data = data %>%
        mutate(label = ifelse(rank(`LFC DMSO/ETP`) <= label_top, `GENE SYMBOL`, NA))
    ggplot(data, aes(x=DMSO, y=`LFC DMSO/ETP`)) +
        geom_point(alpha=0.5) +
        loess_sd(`LFC DMSO/ETP` ~ DMSO, data=data) +
        ggrepel::geom_label_repel(aes(label=label), size=2) +
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

saveRDS(expr, file="overview.rds")
