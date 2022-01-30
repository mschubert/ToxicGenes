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
        viridis::scale_fill_viridis(option="magma") +
        theme_void()
    if (igraph::gsize(net) > 0)
        p = p + geom_edge_link(alpha=0.2)
    p
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'pan/stan-nb.rds'),
        opt('n', 'n_genes', 'number of genes to have non-zero scores', '500'),
        opt('p', 'plotfile', 'pdf', 'bionet.pdf')
    )

    scores = readRDS(args$infile)
    cohort = sub("^([a-zA-Z]+)/.*", "\\1", args$infile)
    ng = as.integer(args$n_genes)

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
    cu = scores %>% filter(z_comp < 0) %>% arrange(z_comp) %>%
        head(ng) %>% top_n(1, z_comp) %>% pull(z_comp) %>% abs() %>% max(1)
    scores2 = -scores$z_comp %>% setNames(scores$gene)
    scores2 = scores2[!is.na(scores2)]
    scores2[scores2 < cu] = 0

    res = as_tbl_graph(runFastHeinz(net, scores2)) %N>%
        left_join(scores %>% dplyr::rename(name=gene))

    pdf(args$plotfile)
    print(plot_net(res, aes(fill=pmax(estimate, -1), size=n_aneup)) +
          ggtitle(sprintf("%s (cutoff: %.2f)", cohort, cu)))
    dev.off()
})
