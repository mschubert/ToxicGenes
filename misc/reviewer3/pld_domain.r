library(dplyr)
library(ggplot2)
gset = import('genesets')
cm = import('../../report/common')

# from: https://www.pnas.org/doi/abs/10.1073/pnas.1109434108, table sd02.xls
tab = readxl::read_xls("sd02.xls")
rrm = rev(rev(tab$`Gene Name`)[-c(1,2)]) # 213x RNA recognition motif
prion = c("DAZ1", "DAZ2", "DAZ3", "DAZAP1", "EWSR1", "FUS", "HNRNPA0",
          "HNRNPA1", "HNRNPA2B1", "RBM14", "TAF15", "TARDBP", "TIA1")

all = cm$get_tox()$`Pan-Cancer` |> pull(gene)
sets = list(
    Compensated = cm$get_comp_genes(pan=TRUE),
    Toxic = cm$get_tox()$`Pan-Cancer` |> filter(is_toxic) |> pull(gene),
    ARGOS = cm$get_argos(pan=TRUE)
)

res = list(
    `RRM over all` = gset$test_fet(valid=all, hits=rrm, sets=sets),
    `PLD over RRM` = gset$test_fet(valid=rrm, hits=prion, sets=sets)
) |> bind_rows(.id="Comparison")

intersect(setdiff(rrm, prion), sets$Compensated)
intersect(prion, sets$Compensated)

saveRDS(res, file="pld_domain.rds")

pdf("pld_domain.pdf", 5, 3)
ggplot(res, aes(x=estimate, y=-log10(p.value))) +
    geom_errorbarh(aes(xmin=conf.low, xmax=conf.high), alpha=0.5) +
    geom_point(aes(shape=Comparison, fill=label), size=2) +
    scale_x_log10() +
    scale_shape_manual(values=c(`RRM over all`=21, `PLD over RRM`=23)) +
    scale_fill_discrete(guide=guide_legend(override.aes=list(shape=21)), name="Gene set") +
    geom_hline(yintercept=-log10(0.05), linetype="dashed") +
    annotate("text", x=0.1, y=-log10(0.05), label="P = 0.05", vjust=-1, hjust=0.8, size=3) +
    labs(x = "Odds ratio (fold enrichment)")
dev.off()
