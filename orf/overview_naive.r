library(ggplot2)
io = import('io')
sys = import('sys')

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

plot_overview_corrected = function(data, title, label_top=20) {
    data = data %>%
        mutate(label = ifelse(rank(z_LFC) <= label_top, `GENE SYMBOL`, NA))
    ggplot(data, aes(x=DMSO, y=z_LFC)) +
        geom_point(alpha=0.5) +
        loess_sd(data$DMSO, data$z_LFC) +
#        stat_loess_sd(color="red", size=1) +
        ggrepel::geom_label_repel(aes(label=label), size=2) +
        labs(title = title,
             x = "DMSO",
             y = "LFC DMSO/ETP z-score") +
        theme_classic()
}

if (is.null(module_name())) {
    args = sys$cmd$parse(
        opt('t', 'tissues', 'txt', '../data/orf/tissues.txt'),
        opt('l', 'orflib', 'txt', '../data/orf/20101003_ORF-size-plasmid.txt'), # ignored (for now)
        opt('s', 'orfscreen', 'xlsx', '../data/orf/ORF_DMSO-ETP_2019-07.xlsx'),
        opt('x', 'orfscreen2', 'txt', '../data/orf/ORF_DMSO_2019-02.txt'),
        opt('o', 'outfile', 'rds', 'overview.rds'),
        opt('p', 'plotfile', 'pdf', 'overview_naive.pdf'))

    tissues = io$read_table(args$tissues, header=TRUE)

    missing = readr::read_tsv(args$orfscreen2) %>%
        tidyr::gather("condition", "value", -(`Construct Barcode`:`BEST GENE MATCH`)) %>%
        filter(condition %in% c("WM266-4 LFC")) # "WM266-4 DMSO" in both

    expr = readxl::read_xlsx(args$orfscreen) %>%
        tidyr::gather("condition", "value", -(`Construct Barcode`:`BEST GENE MATCH`)) %>%
        bind_rows(missing) %>%
        filter(condition != "Hs936T ETP") %>% # don't have LTP/LFC for this line
        mutate(condition = sub("( ORF)?[_-]DMSO", " DMSO", condition),
               condition = sub("LFC$", "LFC DMSO/ETP", condition),
               cells = sub(" .*$", "", condition),
               cells = sub("-ETP", "", cells, fixed=TRUE),
               cells = sub("Kura", "Kuramochi", cells, fixed=TRUE),
               cells = sub("SKBr3", "SK-BR-3", cells, fixed=TRUE),
               cells = sub("LnCAP", "LnCaP", cells, fixed=TRUE),
               cells = sub("SKNEP1", "SK-NEP-1", cells, fixed=TRUE),
               condition = sub("^[A-Za-z0-9-]+ ", "", condition)) %>%
        tidyr::spread(condition, value) %>%
        group_by(cells) %>%
        mutate(z_LFC = loess_z(DMSO, `LFC DMSO/ETP`)) %>%
        ungroup() %>%
        left_join(tissues %>% select(-comment))

    percell = expr %>%
        mutate(cells = sprintf("%s (%s)", cells, tissue)) %>%
        group_by(cells) %>%
        tidyr::nest() %>%
        mutate(plot = purrr::map2(data, cells, plot_overview))

    pdf(args$plotfile)
    for (p in percell$plot)
        print(p)
    dev.off()

    saveRDS(expr, file=args$outfile)
}
