library(dplyr)
sys = import('sys')

avg_undirected = function(genes, ..., resample=1000) {
    avg_ng = function(n) tcga %>% sample_n(n) %>% pull(statistic) %>% mean()

    obs = tcga %>% filter(name %in% genes)
    expected = replicate(resample, avg_ng(length(genes)))
    z = (mean(obs$statistic) - mean(expected)) / sd(expected)
    tibble(observed=mean(obs$statistic), expected=list(expected), z=z,
           p.value=2*pnorm(abs(z), lower.tail=FALSE))
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
    tibble(Gene=obs$name, observed=obs$statistic, expected=list(expected), z=z,
           p.value=2*pnorm(abs(z), lower.tail=FALSE))
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
    tibble(Gene=obs$name, observed=obs$statistic, expected=list(expected), z=z,
           p.value=2*pnorm(z, lower.tail=!as.logical(signs[negpos]+1)))
}

run = function(brc, fun) {
    re = brc %>%
        mutate(res = clustermq::Q(fun, genes=data, negpos=negpos,
            export=list(tcga=tcga), pkgs="dplyr", n_jobs=10)) %>%
        tidyr::unnest("res") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p)
}

plot_top = function(res) {
    res2 = res %>%
        head(10) %>%
        mutate(pos = sprintf("%s:%i-%i", Chr, Peak.Start, Peak.End),
               data = sapply(data, function(s) paste(s, collapse="\n")))
    if ("Gene" %in% colnames(res2))
        res2$data = res2$Gene
    dists = res2 %>%
        select(pos, negpos, expected) %>%
        tidyr::unnest()
    ggplot(res2, aes(x=pos, y=observed)) +
        geom_violin(data=dists, aes(y=expected, fill=negpos)) +
        geom_point(size=3, color="red") +
        geom_label(aes(label=data), hjust=0, size=2, label.size=NA, fill="#ffffffc0") +
        geom_text(aes(label=sprintf("%.2g ", adj.p), hjust=1)) +
        theme(axis.text.x = element_text(angle=30, hjust=1))
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

    pdf("do.pdf", 10, 6)
    print(plot_top(au) + ggtitle("average undirected"))
    print(plot_top(mu) + ggtitle("max undirected"))
    print(plot_top(md) + ggtitle("max directed"))
    dev.off()
})
