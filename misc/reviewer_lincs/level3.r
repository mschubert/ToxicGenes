library(dplyr)
library(ggplot2)
io = import('io')

genes = readr::read_tsv("geneinfo_beta.txt") |> filter(feature_space == "landmark")
meta = readr::read_tsv("instinfo_beta.txt") |>
    filter(pert_type %in% c("ctl_untrt", "trt_oe"), # ctl_vehicle no HEK293T cells
           pert_dose %in% c("200") | is.na(pert_dose), # 300: few genes
           pert_time == "48", # 72/96h: few genes
           qc_pass == 1)

meta1 = meta[meta$pert_type == "ctl_untrt",]
meta2 = meta[meta$pert_type == "trt_oe",]

mat = io$load_gctx("level3_beta_ctl_n188708x12328.gctx")
ctl = mat@mat[as.character(genes$gene_id), meta1$sample_id] |>
    cbind(genes["gene_symbol"]) |>
    tidyr::pivot_longer(-gene_symbol, names_to="sample_id", values_to="expr") |>
    inner_join(meta) |>
    filter(gene_symbol %in% genes$gene_symbol) |>
    select(gene_symbol, cell_iname, pert_itime, expr)

mat = io$load_gctx("level3_beta_trt_oe_n131668x12328.gctx")
dset = mat@mat[as.character(genes$gene_id), meta2$sample_id] |>
    cbind(genes["gene_symbol"]) |>
    tidyr::pivot_longer(-gene_symbol, names_to="sample_id", values_to="expr") |>
    inner_join(meta) |>
    filter(gene_symbol == cmap_name) |>
    select(gene_symbol, cell_iname, pert_idose, pert_itime, expr)

common = inner_join(ctl |> select(-expr) |> distinct(),
                    dset |> filter(!is.na(pert_idose)) |> select(-expr) |> distinct())
long = list(ctl=ctl, oex=dset) |> lapply(inner_join, y=common) |> bind_rows(.id="type")

pdf("level3.pdf", 5, 12)
ggplot(long, aes(x=forcats::fct_reorder(gene_symbol, -expr), y=expr, fill=type)) +
    geom_boxplot(outlier.shape=NA, position=position_dodge(width=0.4), alpha=0.7) +
    facet_wrap(~ cell_iname + pert_itime) +
    coord_flip() +
    labs(x = "gene")
dev.off()
