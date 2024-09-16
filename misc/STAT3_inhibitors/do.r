library(modules)
library(dplyr)
library(ggplot2)
cm = import('../../report/common')

parse_tables = function(sheet, range_index, range_values) {
    idx = readxl::read_xlsx(fname, sheet, range=range_index) |>
        tidyr::fill(...1, ...2, .direction="down") |>
        select(IR=...1, oex=...2, `C188-9_conc`=`C188-9`, HJC052_conc=HJC052) |>
        filter(!is.na(oex))
    val = readxl::read_xlsx(fname, sheet, range=range_values) |>
        filter(!is.na(...1)) |>
        select(`C188-9_value`=`C188-9`, HJC052_value=HJC052) |>
        mutate()
    cbind(idx, val) |>
        tidyr::pivot_longer(c(-IR, -oex), names_to=c("drug", ".value"), names_pattern="(.+)_(.+)")
}

load_sheet = function(sheet) {
    ir0 = parse_tables(sheet, "O41:W50", "A40:M50")
    ir1 = parse_tables(sheet, "AN40:AV49", "Z40:AL50")
    bind_rows(ir0, ir1)
}

fname = "RBM14 data（calculation）Normalize to Luciferase.xlsx"
sheets = readxl::excel_sheets(fname) |> setdiff("Sum up")
dset = sapply(sheets, load_sheet, simplify=FALSE) |>
    bind_rows(.id="sample") |>
    mutate(cell_line = sub("([^ ]+) ([0-9]+)", "\\1", sample),
           rep = sub("([^ ]+) ([0-9]+)", "\\2", sample))

dset2 = dset |>
    group_by(cell_line, drug, oex, IR) |>
    mutate(value = value / value[conc == "0uM"]) |>
    filter(conc == "1uM")

pdf("Rplots.pdf", 6, 3)
ggplot(dset2, aes(x=oex, fill=oex, y=value)) +
    geom_boxplot(alpha=0.7) +
    geom_line(aes(group=paste(sample, drug, IR), linetype=IR), alpha=0.15) +
    ggbeeswarm::geom_quasirandom(aes(shape=drug), width=0.2, alpha=0.9) +
    facet_grid(. ~ cell_line) +
    ggsignif::geom_signif(y_position=1.45, color="black", tip_length=0,
        test=function(...) t.test(..., paired=TRUE),
        map_signif_level=cm$fmt_p, parse=TRUE,
        comparisons=list(c("Luc", "RBM14"))) +
    scale_shape_manual(values=c(21,22)) +
    coord_cartesian(ylim=c(0.5, 1.55)) +
    scale_fill_brewer(palette="Set1") +
    labs(x = "Gene overexpressed",
         shape = "Drug",
         fill = "Gene\noverexpressed") +
    cm$theme_minimal()
dev.off()
