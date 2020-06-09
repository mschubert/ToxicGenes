library(dplyr)
sys = import('sys')

avg_undirected = function(genes, .ignored) {
    avg_ng = function(n) tcga %>% sample_n(n) %>% pull(statistic) %>% mean()

    obs = tcga %>% filter(name %in% genes)
    stat_exp = replicate(1000, avg_ng(length(genes)))
    z = (mean(obs$statistic) - mean(stat_exp)) / sd(stat_exp)
    data.frame(Gene=list(obs$name), z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
}

avg_directed = function(genes, negpos) {
    avg_ng = function(n) tcga %>% sample_n(n) %>% pull(statistic) %>% mean()

    signs = setNames(c(-1,1),c("n","p"))
    obs = tcga %>% filter(name %in% genes) %>% pull(statistic) %>% mean()
    stat_exp = replicate(1000, avg_ng(length(genes)))
    z = (mean(obs$statistic) - mean(stat_exp)) / sd(stat_exp)
    data.frame(Gene=list(obs$name), z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
}

max_undirected = function(genes, .ignored) {
    max_ng = function(n) {
        tcga %>%
            sample_n(n) %>%
            mutate(rnk = rank(-abs(statistic))) %>%
            filter(rnk == 1) %>%
            pull(statistic)
    }

    obs = tcga %>% filter(name %in% genes) %>%
        mutate(rnk = rank(-abs(statistic))) %>%
        filter(rnk == 1)
    stat_exp = replicate(1000, max_ng(length(genes)))
    z = (obs$statistic - mean(stat_exp)) / sd(stat_exp)
    data.frame(Gene=obs$name, z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
}

max_directed = function(genes, negpos) {
    max_ng = function(n) {
        tcga %>%
            sample_n(n) %>%
            mutate(rnk = rank(-signs[negpos]*(statistic))) %>%
            filter(rnk == 1) %>%
            pull(statistic)
    }

    signs = setNames(c(-1,1),c("n","p"))
    obs = tcga %>% filter(name %in% genes) %>%
        mutate(rnk = rank(-signs[negpos]*(statistic))) %>%
        filter(rnk == 1)
    stat_exp = replicate(1000, max_ng(length(genes)))
    z = (obs$statistic - mean(stat_exp)) / sd(stat_exp)
    data.frame(Gene=obs$name, z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
}

run = function(brc, fun) {
    brc %>%
        mutate(res = purrr::map2(data, negpos, fun)) %>%
        tidyr::unnest("res") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p)
}

sys$run({
    tcga = readxl::read_xlsx("../../tcga/pan/rlm3_pur.xlsx", sheet="amp")

    brc = readr::read_tsv("all_BrISCUT_results_50under_200424.txt") %>%
        select(-X1) %>%
        select(type, Chr, Peak.Start, Peak.End, negpos, Gene) %>%
        filter(Gene %in% tcga$name) %>%
        nest_by(type, Chr, Peak.Start, Peak.End, negpos) %>%
        ungroup() %>%
        filter(type == "PANCAN")

    ud = run(brc, max_undirected)
    di = run(brc, max_directed)
})
