library(dplyr)
library(ggplot2)

ov = readRDS("../../model_orf/overview.rds")
ca1 = readRDS("CRISPRa_tox.rds") |> select(`GENE SYMBOL`=gene, ca_HT29=lfc_CRISPRa)
ca2 = readRDS("CRISPRa_Veronica.rds") |> select(`GENE SYMBOL`=gene, ca_BT869=avg_lfc)

wide = ov |>
    group_by(cells, `GENE SYMBOL`) |>
        summarize(`LFC DMSO/ETP` = mean(`LFC DMSO/ETP`)) |>
    ungroup() |>
    select(`GENE SYMBOL`, cells, `LFC DMSO/ETP`) |>
    tidyr::pivot_wider(names_from=cells, values_from=`LFC DMSO/ETP`) |>
    inner_join(ca1) |>
    inner_join(ca2)
mat = data.matrix(wide[-1])
rownames(mat) = wide[[1]]

cmat = reshape2::melt(cor(mat)) |>
    filter(Var1 != Var2) |>
    mutate(type = ifelse(grepl("^ca_", paste(Var1, Var2)), "CRISPRa-ORF", "ORF"))

pdf("crispra_orf_cor.pdf", 5, 4)
ggplot(cmat, aes(x=type, y=value)) +
    ggbeeswarm::geom_quasirandom()
dev.off()
