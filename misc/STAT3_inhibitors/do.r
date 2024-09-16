library(dplyr)
library(ggplot2)

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
