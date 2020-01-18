library(BioNet)
library(ggraph)
library(tidygraph)
library(patchwork)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')

#' Return BioNet Steiner Tree subnetwork
#'
#' @param g  tidygraph-compatible network
#' @param assocs  data.frame with fields: n_smp, p.value, adj.p
#' @param thresh  p-value/fdr cutoff
#' @return  tidygraph object
bionet = function(g, assocs, thresh=0.05, var="adj.p") {
    g = g %>% activate(nodes) %>% left_join(assocs)
    scores = setNames(pull(g, !! rlang::sym(var)), pull(g, name))
    scores[is.na(scores)] = 1
    scores = pmax(-log10(scores) + log10(thresh), 0)
    as_tbl_graph(runFastHeinz(g, scores)) %>%
        activate(edges) %>%
        filter(from != to)
}

#' Plot a network
#'
#' @param net  ggraph-compatible network object
#' @param node_aes  aesthetics mapping for geom_node_point
#' @return  ggplot2 object
plot_net = function(net, node_aes, ...) {
    set.seed(121979) # same layout if same nodes
    p = ggraph(net) +
        geom_node_point(node_aes, ...) +
        geom_node_text(aes(label = name), size=2, repel=TRUE) +
        scale_color_gradient2(low="blue", high="red", mid="black") +
        theme_void()
    if (igraph::gsize(net) > 0)
        p = p + geom_edge_link(alpha=0.2)
    p
}

#' Plot the highest degree nodes of network
#'
#' @param net  ggraph-compatible network object
#' @return  ggplot2 object
plot_col = function(net, var) {
    vs = rlang::sym(var)
    measures = net %N>%
        mutate(degree = centrality_degree()) %>%
        arrange(- !! vs) %>%
        as_tibble() %>%
        head(20) %>%
        filter(!! vs > 0) %>%
        mutate(name = factor(name, levels=rev(name)))

    ggplot(measures, aes_string(x="name", y=var)) +
        geom_col(aes(fill=log10_ks_p)) +
        scale_fill_gradient2(low="blue", high="red", mid="black") +
        coord_flip() +
        guides(fill=FALSE)
}

sys$run({
    args = sys$cmd$parse(
        opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
        opt('f', 'fit', 'fit type id', 'rlm3'),
        opt('b', 'briscut', 'txt file', '../data/all_BrISCUT_results.txt'),
        opt('r', 'ranks', 'xlsx', 'rank_top/pan_rlm3.xlsx'),
        opt('p', 'plotfile', 'pdf', 'bionet/pan/rlm3.pdf'))

    if (args$tissue == "BRCA") {
        fdr = 0.23
    } else if (args$tissue %in% c("pan", "SKCM", "OV")) {
        fdr = 0.25
    } else {
        fdr = 0.3
    }

    read_stats = function(sheet, ranks=args$ranks) {
        tab = readxl::read_xlsx(ranks, sheet=sheet) %>%
            mutate(pseudo_p = 10^(-score))
    }
    top = sapply(readxl::excel_sheets(args$ranks), read_stats, simplify=FALSE)

    cnv = setNames(c(1, -1), c("amp", "del"))
    genomic = readr::read_tsv(args$briscut) %>%
        transmute(cohort=type, name=Gene, log10_ks_p=pmin(log10_ks_p, 10) * cnv[direction])
    if (args$tissue == "pan") {
        genomic = genomic %>% group_by(name) %>% summarize(log10_ks_p = mean(log10_ks_p))
    } else {
        genomic = genomic %>% filter(cohort == args$tissue)
    }

    net = OmnipathR::import_AllInteractions() %>%
        OmnipathR::interaction_graph() %>%
        as_tbl_graph(directed=FALSE) %E>%
        select(-dip_url, -sources, -references) %N>% # ggraph issue #214
        left_join(genomic)

    res = lapply(top, function(t) bionet(net, t, thresh=fdr, var="pseudo_p"))

    pdf(args$plotfile, 8, 7)
    for (i in seq_along(res)) {
        p1 = plot_net(res[[i]] %N>% mutate(degree = centrality_degree()),
                      aes(size=score, alpha=degree, color=log10_ks_p))
        p2 = plot_col(res[[i]], "score")
        p3 = plot_col(res[[i]], "degree")
        print(p1 + ggtitle(names(res)[i]) + (p2 / p3) + plot_layout(widths=c(8,1)))
    }
    dev.off()
})
