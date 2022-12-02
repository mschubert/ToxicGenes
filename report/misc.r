library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
cm = import('./common')

set_cor_fet = function(all) {
    do_fet = function(needle, haystack) {
        both = length(intersect(needle, haystack))
        mat = matrix(c(both,
                       length(setdiff(needle, both)),
                       length(setdiff(haystack, both)),
                       bg - both), nrow=2)
        broom::tidy(fisher.test(mat)) %>% select(estimate, p.value)
    }
    gset = import('genesets')
    go = gset$get_human("CORUM_all")
    set_comp = unique(na.omit(all$gene[all$est_tcga < -0.3 & all$est_ccle < -0.3]))
    set_orf = unique(na.omit(all$gene[all$stat_orf < -5]))
    bg = length(intersect(unlist(go), unlist(all$gene)))

    res = tibble(label = names(go)) %>%
        rowwise() %>%
            mutate(size = length(go[[label]]),
                   comp = do_fet(set_comp, go[[label]]),
                   orf = do_fet(set_orf, go[[label]])) %>%
        ungroup() %>%
        tidyr::unnest_wider(c(comp, orf), names_sep=".")

    plt$denspt(res, aes(x=comp.estimate, y=orf.estimate, label=label, size=size)) +
        scale_size_area(max_size=8, breaks=c(10,100,500,1000), name="Genes in set")
}

corum_cors = function() {
    orf = readxl::read_xlsx("../orf/pan/CORUM_all.xlsx")
    ccle_go = readRDS("../ccle/pan/stan-nb/all/CORUM_all.rds")$amp %>%
        select(label, stat_ccle=statistic)
    tcga_go = readRDS("../tcga/pan/stan-nb_puradj/all/CORUM_all.rds")[[1]] %>%
        select(label, stat_tcga=statistic, size_used)
    both = inner_join(ccle_go, tcga_go) %>% filter(size_used < 1000) %>%
        inner_join(orf %>% select(label=name, stat_orf=statistic)) %>%
        mutate(tcga_ccle = (stat_ccle+stat_tcga)/2)

    m = broom::glance(lm(tcga_ccle ~ stat_orf, data=both))
    lab = sprintf("R^2~`=`~%.2f~\n~p~`=`~%.1g", m$adj.r.squared, m$p.value) %>%
        sub("e", "%*%10^", .)

    plt$denspt(both, aes(x=tcga_ccle, y=stat_orf, label=label), size=size_used) +
        scale_size_area(max_size=16, breaks=c(10,100,500,1000), name="Genes in set") +
        theme_minimal() +
        annotate("text", x=3, y=-8, color="blue", label=lab, parse=TRUE)
}

all = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv")

pdf("misc.pdf", 12, 10)
set_cor_fet(all)
corum_cors()
dev.off()
