library(dplyr)
sys = import('sys')

avg_undirected = function(genes, ..., resample=1000) {
    avg_ng = function(n) tcga %>% sample_n(n) %>% pull(statistic) %>% mean()

    obs = tcga %>% filter(name %in% genes)
    expected = replicate(resample, avg_ng(length(genes)))
    z = (mean(obs$statistic) - mean(expected)) / sd(expected)
    data.frame(z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
}

max_undirected = function(genes, ..., resample=1000) {
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
    expected = replicate(resample, max_ng(length(genes)))
    z = (obs$statistic - mean(expected)) / sd(expected)
    data.frame(Gene=obs$name, z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
}

max_directed = function(genes, negpos, resample=1000) {
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
    expected = replicate(resample, max_ng(length(genes)))
    z = (obs$statistic - mean(expected)) / sd(expected)
    data.frame(Gene=obs$name, z=z, p.value=2*pnorm(abs(z), lower.tail=FALSE))
}

run = function(brc, fun) {
    re = brc %>%
        mutate(res = clustermq::Q(fun, genes=data, negpos=negpos,
                                  export=list(tcga=tcga), n_jobs=10)) %>%
        tidyr::unnest("res") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p)
}

sys$run({
    tcga = readxl::read_xlsx("../../tcga/pan/rlm3_pur.xlsx", sheet="amp") %>%
        filter(!is.na(statistic)) %>%
        select(name, statistic)

    brc = readr::read_tsv("all_BrISCUT_results_50under_200424.txt") %>%
        select(-X1) %>%
        select(type, Chr, Peak.Start, Peak.End, negpos, Gene) %>%
        filter(Gene %in% tcga$name) %>%
        group_by(type, Chr, Peak.Start, Peak.End, negpos) %>%
            tidyr::nest() %>%
        ungroup() %>%
        mutate(data = lapply(data, unlist, use.names=FALSE)) %>% # df->list<chr>
        filter(type == "PANCAN")

    au = run(brc, avg_undirected)
    mu = run(brc, max_undirected)
    md = run(brc, max_directed)
})
