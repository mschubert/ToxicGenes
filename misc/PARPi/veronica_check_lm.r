library(dplyr)
library(ggplot2)

RBM14 <- read.csv(file = 'CCLE_RBM14expCN.csv') %>%
    select(depmap_id, cell_line_display_name, lineage_1,
           expression = Expression.Public.22Q4.RBM14,
           copy_number = Copy.Number.Public.22Q4.RBM14)
PRISM <- read.csv(file = 'PRISMsecondary.csv') %>%
    filter(passed_str_profiling) %>%
    select(depmap_id, broad_id, screen_id, auc, name, moa, target)

result = inner_join(PRISM, RBM14) %>%
    group_by(broad_id, name, moa, lineage_1, expression, copy_number, cell_line_display_name) %>%
        summarize(auc = mean(auc, na.rm=TRUE)) %>%
    group_by(broad_id, name, moa) %>%
        summarize(by_expr = list(broom::tidy(lm(auc ~ lineage_1 + expression))),
                  by_cn = list(broom::tidy(lm(auc ~ lineage_1 + copy_number)))) %>%
    ungroup() %>%
    tidyr::pivot_longer(c(by_expr, by_cn), names_to="comparison") %>%
    tidyr::unnest(value) %>%
    filter(term %in% c("expression", "copy_number")) %>%
    group_by(term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr"),
           type = case_when(
               moa == "PARP inhibitor" ~ "PARPi",
               p.value < 0.05 ~ "p<0.05",
               TRUE ~ "n.s."
           )) %>%
    arrange(adj.p, p.value)

ggplot(result, aes(x=estimate, y=-log10(p.value), color=type)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label=name), max.overlaps=5) +
    facet_wrap(~ term)
