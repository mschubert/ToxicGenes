library(dplyr)
library(ggplot2)
io = import('io')

meta = readr::read_tsv("siginfo_beta.txt") |>
    filter(pert_type == "trt_oe",
           pert_dose %in% c("200"), # 300: few genes
           pert_time == "48", # 72/96h: few genes
           qc_pass == 1)
mat = io$load_gctx("level5_beta_trt_oe_n34171x12328.gctx")
genes = readr::read_tsv("genes.tsv.txt") |>
    filter(`Entrez ID` %in% rownames(mat@mat),
           Type == "landmark")

dset = cbind(genes["Symbol"], mat@mat[as.character(genes$`Entrez ID`),]) |>
    tidyr::pivot_longer(-Symbol, names_to="sig_id", values_to="beta") |>
    inner_join(meta) |>
    filter(Symbol == cmap_name)

ggplot(dset, aes(x=forcats::fct_reorder(cmap_name, -beta), y=beta)) +
    geom_boxplot() +
    facet_wrap(~ cell_mfc_name + pert_idose + pert_itime) +
    coord_flip() +
    labs(x = "gene")
