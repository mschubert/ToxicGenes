library(dplyr)
library(depmap)
library(ExperimentHub)
sys = import('sys')
plt = import('plot')

calc_assocs = function(dset, field) {
    message(as.character(dset), " @ ", field)
    res = get(dset) %>%
        select(depmap_id, gene_name, dependency) %>%
        filter(!is.na(dependency)) %>%
        inner_join(meta) %>%
        mutate(field := !! rlang::sym(field)) %>%
        group_by(gene_name) %>%
            filter(!is.na(dependency), !is.na(field)) %>%
            filter(n_distinct(lineage) > 1) %>%
            summarize(mod = list(broom::tidy(lm(dependency ~ lineage + field)))) %>%
        tidyr::unnest(mod) %>%
        filter(term == "field") %>%
        mutate(p.adj = p.adjust(p.value, method="fdr")) %>%
        arrange(p.adj, p.value)
}

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'depmap.rds'),
    opt('p', 'plotfile', 'pdf', 'depmap.pdf')
)

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

rnai = depmap::depmap_rnai()
crispr_ko = depmap::depmap_crispr()

idx = tidyr::crossing(tibble(dset = c("rnai", "crispr_ko")),
                      tibble(field = c("rna_expression", "log_copy_number"))) %>%
    rowwise() %>%
    mutate(res = list(calc_assocs(dset, field)),
           plot = list(plt$volcano(res) + ggtitle(paste(dset, field))))

pdf(args$plotfile)
for (p in idx$plot)
    print(p)
dev.off()

saveRDS(idx, args$outfile)
