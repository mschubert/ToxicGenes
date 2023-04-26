library(dplyr)
library(depmap)
library(ExperimentHub)
library(patchwork)
sys = import('sys')
plt = import('plot')

calc_assocs = function(dset, field, cond) {
    message(dset, " @ ", field)

    make_fml = function(y, ...) paste(y, "~", paste(list(...), collapse=" + "))
    if (cond == "naive") {
        fml = make_fml("dependency", "lineage", field)
    } else {
        fml = make_fml("dependency", "lineage", sub("RBM14", cond, field), field)
    }

    dsets[[dset]] %>%
        select(depmap_id, group, label, dependency) %>%
        filter(!is.na(dependency)) %>%
        inner_join(meta) %>%
        group_by(group, label) %>%
            filter(n_distinct(lineage) > 1) %>%
            summarize(mod = list(broom::tidy(lm(as.formula(fml))))) %>%
        ungroup() %>%
        tidyr::unnest(mod) %>%
        filter(term == field) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

get_meta = function() {
    tpm = depmap::depmap_TPM() %>%
        filter(gene_name %in% c("RBM14", "CCND1")) %>%
        select(depmap_id, gene_name, rna_expression) %>%
        tidyr::pivot_wider(names_from=gene_name, names_prefix="expr_", values_from=rna_expression)
    copy = depmap::depmap_copyNumber() %>%
        filter(gene_name %in% c("RBM14", "CCND1")) %>%
        select(depmap_id, gene_name, log_copy_number) %>%
        tidyr::pivot_wider(names_from=gene_name, names_prefix="copy_", values_from=log_copy_number)
    depmap::depmap_metadata() %>%
        select(depmap_id, cell_line, lineage, sample_collection_site, primary_or_metastasis, sex) %>%
        inner_join(tpm) %>%
        inner_join(copy)
}

get_dsets = function() {
    prism1 = depmap::depmap_drug_sensitivity() %>%
        mutate(group=compound, label=name) %>%
        split(.$screen_id) %>% setNames(paste0("prism1_", names(.)))
    prism2 = readr::read_csv("PRISMsecondary.csv") %>%
        mutate(group=broad_id, label=name, dependency=auc) %>%
        split(.$screen_id) %>% setNames(paste0("prism2_", names(.)))

    c(list(
        rnai = depmap::depmap_rnai() %>% mutate(group=gene, label=gene_name),
        crispr_ko = depmap::depmap_crispr() %>% mutate(group=gene, label=gene_name)
    ), prism1, prism2)
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'rds', 'depmap.rds'),
        opt('p', 'plotfile', 'pdf', 'depmap.pdf')
    )

    meta = get_meta()
    dsets = get_dsets()

    idx = tidyr::expand_grid(dset = names(dsets),
                             field = c("expr_RBM14", "copy_RBM14"),
                             cond = c("naive", "CCND1")) %>%
        rowwise() %>%
        mutate(res = list(calc_assocs(dset, field, cond)))

    plots = idx %>%
        mutate(plot = list(plt$volcano(res) + ggtitle(sprintf("%s %s (%s)", dset, field, cond)))) %>%
        group_by(dset) %>%
        summarize(asm = list(wrap_plots(plot)))

    pdf(args$plotfile, 12, 12)
    for (p in plots$asm)
        print(p)
    dev.off()

    saveRDS(idx, args$outfile)
})
