import_package("dplyr", attach=TRUE)
import_package("ggplot2", attach=TRUE)
import_package("BioNet", attach=TRUE)
import_package("ggraph", attach=TRUE)
import_package("tidygraph", attach=TRUE)
sys = import('sys')

#' Plot a network
#'
#' @param net  ggraph-compatible network object
#' @param node_aes  aesthetics mapping for geom_node_point
#' @return  ggplot2 object
plot_net = function(net, node_aes, ...) {
    set.seed(121979) # same layout if same nodes
    p = ggraph(net, layout="stress") +
        geom_node_point(node_aes, alpha=0.7, ..., shape=21) +
        geom_node_text(aes(label = name), size=2, repel=TRUE) +
        viridis::scale_fill_viridis(option="magma", direction=-1) +
        theme_void()
    if (igraph::gsize(net) > 0)
        p = p + geom_edge_link(alpha=0.2)
    p
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'overview_24h.rds'),
#        opt('n', 'n_genes', 'number of genes to have non-zero scores', '500'),
        opt('p', 'plotfile', 'pdf', 'bionet_24h.pdf')
    )

    # splice
    scores1 = readRDS(args$infile) %>%
        filter(cline != "all", collection == "genes", stype != "gex", adj.p < 1e-4) %>%
        arrange(adj.p, -abs(statistic)) %>% filter(!duplicated(label))# %>% head(300)
    # diff expr
    scores2 = readRDS(args$infile) %>%
        filter(cline == "all", collection == "genes", junction == "covar", adj.p < 0.3) %>%
        arrange(adj.p) %>% filter(!duplicated(label))

    op = OmnipathR::import_all_interactions()
    net = op %>%
        OmnipathR::interaction_graph() %>%
        as_tbl_graph() %>%
            convert(to_undirected, .clean=TRUE) %>%
            convert(to_simple, .clean=TRUE) %E>%
        select(-.orig_data)
    conn = igraph::components(net)$membership == 1 # doesn't work in filter above
    net = net %N>% filter(conn) %>% activate(edges)

    # calculate cutoff as z>=1 or top <ng> genes
#    cu = scores %>% filter(z_comp < 0) %>% arrange(z_comp) %>%
#        head(ng) %>% top_n(1, z_comp) %>% pull(z_comp) %>% abs() %>% max(1)
#    scores2 = -scores$z_comp %>% setNames(scores$gene)
#    scores2 = scores2[!is.na(scores2)]
#    scores2[scores2 < cu] = 0

    scores1_net = abs(scores1$statistic) %>% setNames(scores1$label)
    scores2_net = abs(scores2$statistic) %>% setNames(scores2$label)

    res1 = as_tbl_graph(runFastHeinz(net, scores1_net)) %N>%
        left_join(scores1 %>% dplyr::rename(name=label))
    res2 = as_tbl_graph(runFastHeinz(net, scores2_net)) %N>%
        left_join(scores2 %>% dplyr::rename(name=label))

    pdf(args$plotfile)
    print(plot_net(res1, aes(fill=score, size=score)) + ggtitle("splice 24h"))
    print(plot_net(res2, aes(fill=score, size=score)) + ggtitle("gex 24h"))
    dev.off()
})
