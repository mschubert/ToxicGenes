library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
sys = import('sys')

args = sys$cmd$parse(
    opt('c', 'ccle', 'xlsx', './ccle/pan_rlm/genes.xlsx'),
    opt('d', 'ccle', 'xlsx', './ccle/pan_rank/genes.xlsx'),
    opt('', '', '', ''))

orf = readxl::read_xlsx("orf/pan/genes.xlsx") %>%
    mutate(fit = "lm")

ccle = list(
    rlm = readxl::read_xlsx("ccle/pan_rlm/genes.xlsx"),
    rank = readxl::read_xlsx("ccle/pan_rank/genes.xlsx")
) %>% bind_rows(.id="fit")

tcga_naive = list(
    rlm = readxl::read_xlsx("./tcga/naive/pan_rlm/genes.xlsx"),
    rank = readxl::read_xlsx("./tcga/naive/pan_rank/genes.xlsx")
) %>% bind_rows(.id="fit")
tcga_pur = list(
    rlm = readxl::read_xlsx("./tcga/pur/pan_rlm/genes.xlsx"),
    rank = readxl::read_xlsx("./tcga/pur/pan_rank/genes.xlsx")
) %>% bind_rows(.id="fit")
tcga_puradj = list(
    rlm = readxl::read_xlsx("./tcga/puradj/pan_rlm/genes.xlsx"),
    rank = readxl::read_xlsx("./tcga/puradj/pan_rank/genes.xlsx")
) %>% bind_rows(.id="fit")

tcga = list(
    naive = tcga_naive,
    pur = tcga_pur,
    puradj = tcga_puradj
) %>% bind_rows(.id="adj")

dset = list(orf=orf, ccle=ccle, tcga=tcga) %>%
    bind_rows(.id="dset") %>%
    select(name, dset, fit, adj, statistic) %>%
    mutate(adj = ifelse(is.na(adj), "none", adj))

top = dset %>%
    group_by(dset, fit, adj) %>%
    mutate(score = rank(-statistic) / length(statistic)) %>%
    group_by(name, dset, fit) %>%
    summarize(score = mean(score, na.rm=TRUE)) %>%
    group_by(name, dset) %>%
    summarize(score = mean(score, na.rm=TRUE)) %>%
    group_by(name) %>%
    filter(all(score > 0)) %>%
    summarize(n_dset = length(name),
              score = sum(score, na.rm=TRUE) / n_dset^0) %>% # ^0 is sum (CDKN1A), ^1 is mean (LOCxxx)
    arrange(-score)

do_plot = function(gene) {
    dset %>%
        filter(name == gene) %>%
        mutate(dset = relevel(factor(dset), "orf")) %>%
        ggplot(aes(x=name, y = statistic, color=adj)) +
            geom_point(size=5) +
            facet_wrap(~ dset + fit, scale="free_x", nrow=1) +
            expand_limits(y=0)
}

pdf("candidates.pdf")
for (g in head(top$name, 10))
    print(do_plot(g) + ggtitle(g))
dev.off()
