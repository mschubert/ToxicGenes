library(dplyr)
library(ggplot2)
util = import('./overview_naive')

plot_one = function(df, x, y, nbin=50) {
    rx = diff(range(df[[as.character(rlang::ensym(x))]], na.rm=TRUE))
    ry = diff(range(df[[as.character(rlang::ensym(y))]], na.rm=TRUE))
    ggplot(df, aes(x={{ x }}, y={{ y }})) +
        stat_bin2d(aes(fill=after_stat(count)), binwidth=c(rx/nbin, ry/nbin)) +
        scale_fill_distiller(palette="Spectral") +
        geom_smooth(color="red") +
        facet_wrap(~ cline)
}

size = readr::read_tsv("../data/orf/20101003_ORF-size-plasmid.txt")
info = size[c("X9", "X10")]
info = info[rowSums(is.na(info)) < ncol(info),]
size = size %>% select(-WARNINGS, -(X8:X10))

etp = readxl::read_xlsx("../data/orf/ORF_DMSO-ETP_2019-07.xlsx") %>%
    tidyr::gather("condition", "value", -(1:4)) %>%
    mutate(cline = sub("(ETP|(LFC )?DMSO(/ETP)?)", "", condition),
           cline = sub("(LFC|ORF_?)", "", cline),
           cline = sub("-$", "", cline),
           cline = stringr::str_trim(cline),
           cline = sub("LnCAP", "LnCaP", cline, fixed=TRUE),
           cline = sub("SKNEP1", "SK-NEP-1", cline, fixed=TRUE),
           condition = case_when(
               grepl("LFC", condition) ~ "LFC DMSO/ETP",
               grepl("ETP", condition) ~ "early",
               grepl("DMSO", condition) ~ "late"
           )) %>%
    tidyr::spread("condition", "value")

ov = readRDS("overview.rds") %>%
    select(`Construct IDs`, cline=cells, z_LFC)

lib = size %>%
    tidyr::gather("cline", "log rpm", -(`Construct Barcode`:`ORF LENGTH`))

both = inner_join(etp, size) %>%
    left_join(ov)
#plot(both$`INSERT LENGTH`, both$`ORF LENGTH`)

pdf("early-vs-late.pdf", 10, 8)
print(plot_one(lib, `ORF LENGTH`, `log rpm`) + ggtitle("Library representation"))
print(plot_one(both, `ORF LENGTH`, early) + labs(y="log rpm", title="Early time point"))
print(plot_one(both, `ORF LENGTH`, late) + labs(y="log rpm", title="Late time point"))
print(plot_one(both, `ORF LENGTH`, `LFC DMSO/ETP`) + ggtitle("Early vs. Late TP"))
print(plot_one(both, `ORF LENGTH`, z_LFC) + ggtitle("Early vs. Late TP"))
dev.off()
