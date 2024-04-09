library(dplyr)
library(ggplot2)
seq = import('seq')
cm = import('../../report/common')

comp = cm$get_comp_genes(pan=TRUE)

chrs = seq$gene_table() |>
    select(gene=external_gene_name, chr=chromosome_name, loc=start_position) |>
    group_by(gene, chr) |>
        summarize(loc = mean(loc, na.rm=TRUE)) |>
    ungroup()

dist = cm$get_tox()$`Pan-Cancer` |>
    inner_join(chrs) |>
    group_by(chr) |>
        mutate(is_comp = gene %in% comp,
               dist = sapply(loc, function(l) min(abs(l - loc[is_comp])) / 1e3)) |>
    ungroup() |>
    filter(dist < 500)

mods = dist |>
    filter(!is_comp) |>
    summarize(mod = list(broom::tidy(lm(estimate ~ dist)))) |>
    tidyr::unnest(mod) |>
    filter(term == "dist") |>
    mutate(lab = sprintf("P = %.2g", p.value))

pdf("tox_distance.pdf", 5.5, 3.5)
ggplot(dist, aes(x=dist, y=estimate)) +
    geom_point(aes(color=is_comp), alpha=0.5) +
    scale_color_manual(values=c("FALSE"="black", "TRUE"="#de493d"), name="Compensated") +
    geom_text(data=mods, aes(x=50, y=-1, label=lab), hjust=0, color="blue") +
    geom_smooth(data=dist[!dist$is_comp,], method="lm") +
    labs(x = "Distance to closest compensated gene (Kb)",
         y = "log2 fold change in ORF")
dev.off()
