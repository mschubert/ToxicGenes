library(dplyr)
library(ggplot2)
library(DESeq2)
cm = import('../../report/common')

set1 = readr::read_tsv("GCF6748_Calabrese_setA_in_HT29.txt")
colnames(set1) = sub("^a_", "", colnames(set1))
set2 = readr::read_tsv("GCF6748_Calabrese_setB_in_HT29.txt")
colnames(set2) = sub("^b_", "", colnames(set2))

dset = bind_rows(set1, set2)
meta = data.frame(sample = colnames(dset)[-(1:3)]) |>
    mutate(sample = factor(sample) |> relevel("t0_1"),
           treated = as.integer(sample != "t0_1"))
res = DESeqDataSetFromMatrix(dset[-(1:3)], meta, design=~treated) |>
    estimateSizeFactors() |> # size factors on control guides?
    DESeq() |>
    results(name="treated")

res2 = as_tibble(cbind(dset[1:3], as.data.frame(res))) |>
    arrange(padj) |>
    group_by(gene_symbol) |>
    summarize(stat_CRISPRa = mean(stat, na.rm=TRUE)) |>
    dplyr::rename(gene = gene_symbol)

argos = cm$get_argos(pan=TRUE)
tox = cm$get_tox()$`Pan-Cancer` |>
    inner_join(res2) |>
    mutate(label = ifelse(gene %in% argos, gene, NA_character_))

m = broom::tidy(lm(stat_CRISPRa ~ statistic, data=tox))


pdf("CRISPRa_tox.pdf")
ggplot(tox, aes(x=statistic, y=stat_CRISPRa, color=is_toxic)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label=label), color="black") +
    geom_smooth(method="lm", color="blue") +
    ggtitle("NKI HT29") +
    annotate("text", x=-12, y=-4, color="blue",
             label=sprintf("P = %.2g", m$p.value[m$term == "statistic"]))
dev.off()
