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
        mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
        arrange(p.adj, p.value)
}

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'depmap.rds'),
    opt('p', 'plotfile', 'pdf', 'depmap.pdf')
)

tpm = depmap::depmap_TPM() %>%
    filter(gene_name %in% c("RBM14", "CCND1")) %>%
    select(depmap_id, gene_name, rna_expression) %>%
    tidyr::pivot_wider(names_from=gene_name, names_prefix="expr_", values_from=rna_expression)
copy = depmap::depmap_copyNumber() %>%
    filter(gene_name %in% c("RBM14", "CCND1")) %>%
    select(depmap_id, gene_name, log_copy_number) %>%
    tidyr::pivot_wider(names_from=gene_name, names_prefix="copy_", values_from=log_copy_number)
meta = depmap::depmap_metadata() %>%
    select(depmap_id, cell_line, lineage, sample_collection_site, primary_or_metastasis, sex) %>%
    inner_join(tpm) %>%
    inner_join(copy)

drug = depmap::drug_sensitivity_21Q2() %>% mutate(group=compound, label=name)
dsets = list(
    rnai = depmap::depmap_rnai() %>% mutate(group=gene, label=gene_name),
    crispr_ko = depmap::depmap_crispr() %>% mutate(group=gene, label=gene_name),
    drug_hts = drug %>% filter(screen_id == "HTS"),
    drug_mts004 = drug %>% filter(screen_id == "MTS004")
)

idx = tidyr::expand_grid(dset = c("rnai", "crispr_ko", "drug_hts", "drug_mts004"),
                         field = c("expr_RBM14", "copy_RBM14"),
                         cond = c("naive", "CCND1")) %>%
    rowwise() %>%
    mutate(res = list(calc_assocs(dset, field, cond)))

plots = idx %>%
    mutate(plot = list(plt$volcano(res) + ggtitle(sprintf("%s %s (%s)", dset, field, cond)))) %>%
    select(-res) %>%
    tidyr::pivot_wider(names_from="cond", values_from="plot") %>%
    rowwise() %>%
    mutate(plot = list(wrap_plots(naive, CCND1)))

pdf(args$plotfile, 12, 6)
for (p in plots$plot)
    print(p)
dev.off()

saveRDS(idx, args$outfile)
