library(dplyr)
library(patchwork)
library(ggrepel)
sys = import('sys')
plt = import('plot')

plot_one = function(field, cond, res, genes) {
    cur = inner_join(res %>% select(gene=label, dependency=statistic, fdr_dep=adj.p),
                     genes %>% select(gene=gene_name, diff_expr=stat, fdr_de=padj)) %>%
        mutate(group = ifelse(dependency < 0 & fdr_dep < 0.1 & fdr_de < 0.1, "de+dep", "other"))
    m = broom::tidy(lm(dependency ~ diff_expr, data=cur)) %>% filter(term == "diff_expr")
    plt$denspt(cur, aes(x=diff_expr, y=dependency, label=gene, color=group, alpha=0.8),
               draw_pt=1000, draw_label=100) +
        labs(title = sprintf("%s (%s)", field, cond),
             subtitle = sprintf("p=%.2g", m$p.value)) +
        scale_color_manual(values=c("de+dep"="firebrick", other="black"))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'depmap.rds'),
    opt('p', 'plotfile', 'pdf', 'depmap_DEcor.pdf')
)

de_genes = tibble(time = c("8h", "24h", "all")) %>% rowwise() %>%
    mutate(dset = list(readRDS(sprintf("../../data/rnaseq/diff_expr_%s.rds", time)))) %>%
    tidyr::unnest(dset) %>%
    filter(cond == "all_covar") %>%
    select(time, genes)

deps = readRDS(args$infile) %>%
    filter(dset %in% c("rnai", "crispr_ko"), field != "expr_TMEM59L")

dset = tidyr::crossing(deps, de_genes) %>%
    rowwise() %>%
    mutate(plot = list(plot_one(field, cond, res, genes))) %>%
    group_by(time, dset) %>%
    summarize(asm = list((plt$text(paste(dset[1], time[1]), size=8) / wrap_plots(plot)) +
                         plot_layout(heights=c(1,20), guides="collect")))

pdf(args$plotfile, 18, 16)
for (p in dset$asm)
    print(p)
dev.off()
