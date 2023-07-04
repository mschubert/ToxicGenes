library(dplyr)
library(patchwork)
sys = import('sys')
plt = import('plot')
dm = import('./depmap')

get_DEgenes = function() {
    tibble(time = c("8h", "24h", "all")) %>% rowwise() %>%
        mutate(dset = list(readRDS(sprintf("../../data/rnaseq/diff_expr_%s.rds", time)))) %>%
        tidyr::unnest(dset) %>%
        filter(cond %in% c("H838", "H1650", "HCC70", "ZR751")) %>%
        select(time, genes) %>%
        tidyr::unnest(genes) %>%
        group_by(time, gene_name) %>%
        filter(sum(!is.na(stat)) == 4) %>%
        summarize(res = list(broom::tidy(lm(stat ~ 1)))) %>%
        tidyr::unnest(res) %>%
        slice_min(p.value, n=100) %>%
        ungroup()
}

get_scores = function(de_genes) {
    de_genes = get_DEgenes() #%>%
#        mutate(genes = lapply(genes, . %>% slice_max(abs(stat), n=100))) %>%
#        tidyr::unnest(genes)
    scores = depmap::depmap_TPM() %>%
        inner_join(de_genes, relationship="many-to-many") %>%
        group_by(depmap_id, time) %>%
            summarize(score = sum(rna_expression * statistic)) %>%
        group_by(time) %>%
            mutate(score = scale(score)) %>%
        ungroup() %>%
        tidyr::pivot_wider(names_from=time, names_prefix="sig_", values_from=score)

    depmap::depmap_metadata() %>%
        select(depmap_id, cell_line, lineage, sample_collection_site, primary_or_metastasis, sex) %>%
        inner_join(scores)
}

calc_assocs = function(dset, field) {
    message(dset, " @ ", field)
    dsets[[dset]] %>%
        select(depmap_id, group, label, dependency) %>%
        filter(!is.na(dependency)) %>%
        inner_join(scores) %>%
        mutate(field := !! rlang::sym(field)) %>%
        group_by(group, label) %>%
            filter(n_distinct(lineage) > 1) %>%
            summarize(mod = list(broom::tidy(lm(dependency ~ lineage + field)))) %>%
        ungroup() %>%
        tidyr::unnest(mod) %>%
        filter(term == "field") %>%
        select(-term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'depmap.rds'), # ignored, easier snakemake
    opt('p', 'plotfile', 'pdf', 'depmap_DEsig2.pdf')
)

dsets = dm$get_dsets()
scores = get_scores()
idx = tidyr::expand_grid(dset = names(dsets),
                         field = paste0("sig_", c("8h", "24h", "all"))) %>%
    rowwise() %>%
    mutate(res = list(calc_assocs(dset, field)))

plots = idx %>%
    mutate(plot = list(plt$volcano(res) + ggtitle(sprintf("%s %s", dset, field)))) %>%
    mutate(dset = factor(dset, levels=unique(dset))) %>%
    group_by(dset) %>%
    summarize(asm = list(wrap_plots(plot, nrow=1)))

pdf(args$plotfile, 16, 6)
for (p in plots$asm)
    print(p)
dev.off()
