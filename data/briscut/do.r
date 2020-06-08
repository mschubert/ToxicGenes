library(dplyr)

avg_one = function(df) {
    stat = tcga %>% filter(name %in% df$Gene) %>% pull(statistic)
    tt = t.test(stat) %>%
        broom::tidy() %>%
        select(estimate, p.value)
}

max_one = function(genes) {
    max_ng = function(n) {
        tcga %>%
            sample_n(n) %>%
            mutate(rnk = rank(-abs(statistic))) %>%
            filter(rnk == 1) %>%
            pull(statistic)
    }

    stat_obs = tcga %>% filter(name %in% genes) %>%
        mutate(rnk = rank(-abs(statistic))) %>%
        filter(rnk == 1) %>%
        pull(statistic)
    stat_exp = replicate(100, max_ng(length(genes)))
    z = (stat_obs - mean(stat_exp)) / sd(stat_exp)
    #list(z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
    2*pnorm(abs(z), lower.tail=FALSE)
}

tcga = readxl::read_xlsx("../../tcga/pan/rlm3_pur.xlsx", sheet="amp")

brc = readr::read_tsv("all_BrISCUT_results_50under_200424.txt") %>%
    select(-X1) %>%
    select(type, Chr, Peak.Start, Peak.End, negpos, Gene) %>%
    filter(Gene %in% tcga$name) %>%
    nest_by(type, Chr, Peak.Start, Peak.End, negpos) %>%
    filter(type == "PANCAN") %>%
    mutate(maxdevp = purrr::map_dbl(data, max_one)) %>%
    arrange(maxdevp)


