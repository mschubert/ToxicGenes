library(dplyr)
library(ggplot2)
seq = import('seq')
cm = import('../../report/common')

plot_one = function(dists) {
    mods = dists |>
        group_by(Sample) |>
        summarize(mod = list(broom::tidy(lm(LFC ~ dist)))) |>
        tidyr::unnest(mod) |>
        filter(term == "dist") |>
        mutate(lab = sprintf("P = %.2g", p.value))

    ggplot(dists, aes(x=dist, y=LFC)) +
        geom_point(alpha=0.5) +
        geom_text(data=mods, aes(x=0, y=10, label=lab), hjust=0, color="blue") +
        facet_wrap(~ Sample, scales="free_x") +
        geom_smooth(method="lm") +
        labs(x = "Distance to closest compensated gene (Mb)",
             y = "log2 fold change over expected")
}

# get Fig S2e data, calc dist to comp gene, make plot
lookup = c(SS6="chr +7", SS51="+7 +22", SS111="+8 +9 +18")
chrs = seq$gene_table() |>
    select(label=external_gene_name, chr=chromosome_name, loc=start_position) |>
    group_by(label, chr) |>
        summarize(loc = mean(loc, na.rm=TRUE)) |>
    ungroup()
comp = cm$get_comp_genes(pan=TRUE)
dset = readRDS("../../data/rnaseq_rpe1_broad/compute_fcs.rds") |>
    tidyr::unnest(genes) |>
    inner_join(chrs) |>
    filter((term == "SS6" & chr == "7") |
           (term == "SS51" & chr %in% c("7", "22")) |
           (term == "SS111" & chr %in% c("8", "9", "18"))) |>
    transmute(Sample = factor(lookup[term], levels=lookup),
              Gene=label, chr=chr, LFC=log2FoldChange, loc=loc,
              status = ifelse(Gene %in% comp, "Compensated", "Background"),
              status = factor(status, levels=c("Background", "Compensated"))) |>
    group_by(Sample, chr) |>
        mutate(LFC = scale(LFC, scale=FALSE)[,1]) |>
    ungroup()

dists = dset |>
    group_by(Sample, chr) |>
    mutate(dist = sapply(loc, function(l) min(abs(l - loc[status == "Compensated"])) / 1e6))

pdf("gene_distance.pdf", 6, 4)
plot_one(dists)
plot_one(dists |> filter(dist <= 1))
dev.off()
