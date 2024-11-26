library(dplyr)

wt = function(x, file) write.table(x, file, quote=FALSE, row.names=FALSE, col.names=FALSE)

ccle = readxl::read_xlsx("./SuppData1_CCLE-comp.xlsx", sheet="Pan-Cancer")
tcga = readxl::read_xlsx("./SuppData2_TCGA-comp.xlsx", sheet="Pan-Cancer")
orf = readxl::read_xlsx("./SuppData4_ORF-toxicity.xlsx", sheet="Pan-Cancer") |>
    filter(is_toxic) |> pull(gene)

ccle_comp = ccle |> filter(type == "Compensated") |> pull(gene)
tcga_comp = tcga |> filter(type == "Compensated") |> pull(gene)
both = intersect(ccle_comp, tcga_comp)

wt(ccle_comp, "./SuppInfo/SuppNote-CCLEcomp.txt")
wt(tcga_comp, "./SuppInfo/SuppNote-TCGAcomp.txt")
wt(both, "./SuppInfo/SuppNote-comp.txt")
wt(orf, "./SuppInfo/SuppNote-toxic.txt")
wt(intersect(both, orf), "./SuppInfo/SuppNote-ARGOS.txt")
