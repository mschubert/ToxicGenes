library(dplyr)
library(ggplot2)
cm = import('../../report/common')

# CHECK cor of (CRISPRa or ORF) to COMP -> see if they are comparable


comp = cm$get_comp_tissue() |> filter(tissue == "Pan-Cancer") |> select(gene, compensation)
ov = readRDS("../../model_orf/overview.rds")
ca1 = readRDS("CRISPRa_tox.rds") |> select(gene, ca_HT29=lfc_CRISPRa)
ca2 = readRDS("CRISPRa_Veronica.rds") |> select(gene, ca_BT869=avg_lfc)

wide = ov |>
    dplyr::rename(gene = `GENE SYMBOL`) |>
    group_by(cells, gene) |>
        summarize(`LFC DMSO/ETP` = mean(`LFC DMSO/ETP`)) |>
    ungroup() |>
    select(gene, cells, `LFC DMSO/ETP`) |>
    tidyr::pivot_wider(names_from=cells, values_from=`LFC DMSO/ETP`) |>
    inner_join(ca1) |>
    inner_join(ca2) |>
    inner_join(comp) |>
    na.omit()
mat = data.matrix(wide[-1])
rownames(mat) = wide[[1]]

vars = colnames(wide) |> setdiff(c("gene", "compensation"))
cors = lapply(vars, \(v) broom::tidy(cor.test(wide[[v]], wide[["compensation"]]))) |>
    setNames(vars) |>
    bind_rows(.id="sample") |>
    mutate(type = ifelse(grepl("^ca_", sample), "CRISPRa", "ORF")) |>
    mutate(sample = sub("ca_", "", sample))

pdf("crispra_orf_comp_cor.pdf", 5, 4)
ggplot(cors, aes(x=type, y=p.value, fill=type)) +
    geom_boxplot() +
    geom_point() +
    scale_y_log10() +
    ggrepel::geom_text_repel(aes(label=sample), max.overlaps=Inf)
dev.off()
