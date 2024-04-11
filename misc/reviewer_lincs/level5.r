library(dplyr)
library(ggplot2)
io = import('io')

meta = readr::read_tsv("siginfo_beta.txt") |>
    filter(pert_type == "trt_oe",
           pert_dose %in% c("200"), # 300: few genes
           pert_time == "48", # 72/96h: few genes
           qc_pass == 1)
mat = io$load_gctx("level5_beta_trt_oe_n34171x12328.gctx")
genes = readr::read_tsv("geneinfo_beta.txt") |> filter(feature_space == "landmark")

dset = mat@mat[as.character(genes$gene_id), meta$sig_id] |>
    cbind(genes["gene_symbol"]) |>
    tidyr::pivot_longer(-gene_symbol, names_to="sig_id", values_to="zscore") |>
    inner_join(meta) |>
    filter(gene_symbol == cmap_name)

pdf("level5.pdf")
ggplot(dset, aes(x=forcats::fct_reorder(cmap_name, -zscore), y=zscore)) +
    geom_boxplot() +
    facet_wrap(~ cell_mfc_name + pert_idose + pert_itime) +
    coord_flip() +
    labs(x = "gene")
dev.off()
