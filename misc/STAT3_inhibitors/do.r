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

resp_curves = function(dset) {
    ggplot(dset, aes(x=conc, y=value, color=oex, group=paste(sample, oex, IR, rep))) +
        geom_line(aes(linetype=IR)) +
        facet_grid(drug ~ cell_line)
}

luc_vs_rbm_1uM = function(dset) {
    ggplot(dset |> filter(conc == "1uM"), aes(x=oex, fill=oex, y=value)) +
        geom_hline(yintercept=1, linetype="dashed", color="black") +
        geom_boxplot(alpha=0.7) +
        geom_line(aes(group=paste(sample, drug, IR), linetype=IR), alpha=0.15) +
        ggbeeswarm::geom_quasirandom(aes(shape=IR), width=0.2, alpha=0.9) +
        facet_grid(drug ~ cell_line) +
        ggsignif::geom_signif(y_position=1.2, color="black", tip_length=0,
            test=function(...) t.test(..., paired=TRUE),
            map_signif_level=cm$fmt_p, parse=TRUE,
            comparisons=list(c("Luc", "RBM14"))) +
        scale_shape_manual(values=c(21,22)) +
        coord_cartesian(ylim=c(0.5, 1.3)) +
        scale_fill_brewer(palette="Set1") +
        labs(x = "Gene overexpressed",
             y = "Relative viability 1 uM drug vs. DMSO",
             shape = "IR",
             fill = "Gene\noverexpressed") +
        cm$theme_minimal()
}

fname = "RBM14 data（calculation）Normalize to Luciferase.xlsx"
sheets = readxl::excel_sheets(fname) |> setdiff("Sum up")
dset = sapply(sheets, load_sheet, simplify=FALSE) |>
    bind_rows(.id="sample") |>
    mutate(cell_line = sub("([^ ]+) ([0-9]+)", "\\1", sample),
           conc = factor(conc, levels=c("0uM", "0.1uM", "1uM", "10uM")),
           rep = sub("([^ ]+) ([0-9]+)", "\\2", sample)) |>
    group_by(sample, drug, IR) |>
        mutate(value = value / mean(value[conc == "0uM" & oex == "Luc"])) |>
    ungroup()

dset2 = dset |>
    group_by(sample, drug, oex, IR) |>
        mutate(value = value / mean(value[conc == "0uM"])) |>
    ungroup()

# 0 uM in HCC70 is unreliable, plate effect?
dset3 = dset |>
    filter(! (cell_line == "HCC70" & conc == "0uM")) |>
    mutate(conc = case_when(
        cell_line == "HCC70" & conc == "0.1uM" ~ factor("0uM"),
        TRUE ~ conc
    )) |>
    group_by(sample, drug, oex, IR) |>
        mutate(value = value / mean(value[conc == "0uM"])) |>
    ungroup()

pdf("Rplots.pdf", 7, 5)
print(resp_curves(dset) + ggtitle("Viability normalized to sample + Luc"))
print(resp_curves(dset2) + ggtitle("Viability normalized to sample + treatment"))
print(luc_vs_rbm_1uM(dset2))
print(resp_curves(dset3) + ggtitle("Viability normalized to sample + treatment + HCC70 0.1uM"))
print(luc_vs_rbm_1uM(dset3))
dev.off()
