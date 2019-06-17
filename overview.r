library(cowplot)
io = import('io')

tissues = io$read_table("tissues.txt", header=TRUE)

loess_sd = function(x, y) {
    mod = msir::loess.sd(x, y)
    df = data.frame(x = mod$x, y=mod$y, sd=mod$sd)
    df = df[seq(1, nrow(df), length.out=100),]
    list(geom_line(data=df, aes(x=x, y=y), color="red", size=1),
         geom_line(data=df, aes(x=x, y=y+sd), color="red", linetype="dotted", size=1),
         geom_line(data=df, aes(x=x, y=y-sd), color="red", linetype="dotted", size=1))
}

loess_z = function(x, y) {
    mod = msir::loess.sd(x, y)
    df = data.frame(x = mod$x, y=mod$y, sd=mod$sd)
    df = df[rank(x),]
    stopifnot(x == df$x)
    (y - df$y) / df$sd
}

stat_loess_sd = function(mapping = NULL, data = NULL, geom = "line",
                        position = "identity", na.rm = FALSE, show.legend = NA, 
                        inherit.aes = TRUE, ...) {

    statLoessSd = ggproto("statLoessSd", Stat,
        compute_group = function(data, scales) {
            mod = msir::loess.sd(data$x, data$y)
            df = data.frame(x = mod$x, y=mod$y, sd=mod$sd)
            df = df[seq(1, nrow(df), length.out=100),]
        },
        required_aes = c("x", "y")
    )

    layer(stat = statLoessSd, data = data, mapping = mapping, geom = geom, 
          position = position, show.legend = show.legend, inherit.aes = inherit.aes,
          params = list(na.rm = na.rm, ...))
}

plot_overview = function(data, title, label_top=20) {
    data = data %>%
        mutate(label = ifelse(rank(`LFC DMSO/ETP`) <= label_top, `GENE SYMBOL`, NA))
    ggplot(data, aes(x=DMSO, y=`LFC DMSO/ETP`)) +
        geom_point(alpha=0.5) +
        loess_sd(data$DMSO, data$`LFC DMSO/ETP`) +
#        stat_loess_sd(color="red", size=1) +
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
    group_by(cells) %>%
    mutate(z_LFC = loess_z(DMSO, `LFC DMSO/ETP`)) %>%
    ungroup() %>%
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
