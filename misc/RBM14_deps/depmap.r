library(dplyr)
library(depmap)
library(ExperimentHub)
sys = import('sys')

calc_assocs = function() {
    res = dset %>%
        select(depmap_id, gene_name, dependency) %>%
        filter(!is.na(dependency)) %>%
        inner_join(meta) %>%
        group_by(gene_name) %>%
        filter(n_distinct(lineage) > 1) %>%
        summarize(mod = list(broom::tidy(lm(dependency ~ lineage + rna_expression)))) %>%
        tidyr::unnest()
}

tpm = depmap::depmap_TPM() %>%
    filter(gene_name == "RBM14") %>%
    select(depmap_id, rna_expression)
copy = depmap::depmap_copyNumber() %>%
    filter(gene_name == "RBM14") %>%
    select(depmap_id, log_copy_number)
meta = depmap::depmap_metadata() %>%
    select(depmap_id, cell_line, lineage, sample_collection_site, primary_or_metastasis, sex) %>%
    inner_join(tpm) %>%
    inner_join(copy)

screen = depmap::depmap_rnai()
ko = depmap::depmap_crispr()
