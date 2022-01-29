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
    p = ggraph(net) +
        geom_node_point(node_aes, alpha=0.7, ..., shape=21) +
        geom_node_text(aes(label = name), size=2, repel=TRUE) +
        viridis::scale_fill_viridis(option="magma") +
        theme_void()
    if (igraph::gsize(net) > 0)
        p = p + geom_edge_link(alpha=0.2)
    p
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'pan/stan-nb.rds'),
        opt('c', 'cutoff', 'min z_score to include', '4'),
        opt('p', 'plotfile', 'pdf', 'bionet.pdf')
    )

    scores = readRDS(args$infile)
    cu = as.numeric(args$cutoff)

    op = OmnipathR::import_all_interactions()
    net = op %>%
        OmnipathR::interaction_graph() %>%
        as_tbl_graph() %>%
            convert(to_undirected, .clean=TRUE) %>%
            convert(to_simple, .clean=TRUE) %E>%
        select(-.orig_data)
    conn = igraph::components(net)$membership == 1 # doesn't work in filter above
    net = net %N>% filter(conn) %>% activate(edges)

    scores2 = abs(scores$z_comp %>% setNames(scores$gene) %>% pmin(-cu) + cu)
    scores2 = scores2[!is.na(scores2)]
    res = as_tbl_graph(runFastHeinz(net, scores2)) %N>%
        left_join(scores %>% dplyr::rename(name=gene))

    pdf(args$plotfile)
    print(plot_net(res, aes(fill=estimate, size=n_aneup)))
    dev.off()
})
