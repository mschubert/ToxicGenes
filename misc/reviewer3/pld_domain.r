library(dplyr)
library(ggplot2)
idmap = import('process/idmap')
cm = import('../../report/common')

# from: https://www.nature.com/articles/s41598-017-09714-z
pld = readr::read_file("pld_genes.txt")
ensg = stringr::str_extract_all(pld, "ENSG[0-9]+")[[1]]
hg = sort(unique(idmap$gene(ensg, to="hgnc_symbol")))

all = cm$get_tox()$`Pan-Cancer` |> pull(gene)
sets = list(
    Compensated = cm$get_comp_genes(pan=TRUE),
    Toxic = cm$get_tox()$`Pan-Cancer` |> filter(is_toxic) |> pull(gene),
    ARGOS = cm$get_argos(pan=TRUE)
)

res = lapply(sets, function(x) {
    fisher.test(matrix(c(length(intersect(hg, x)), length(setdiff(hg, x)),
                         length(intersect(all, hg)), length(setdiff(all, hg))), nrow=2, ncol=2)) |>
        broom::tidy()
}) |>
    bind_rows(.id="set")
