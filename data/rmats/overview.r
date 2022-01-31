library(dplyr)
library(ggplot2)
sys = import('sys')

sys$run({
    args = sys$cmd$parse(
        opt('d', 'diff_expr', 'rds', 'diff_expr_24h.rds'),
        opt('s', 'splice_all', 'rds', './splice/paired-all_rbm24_vs_luc24.rds'),
        opt('p', 'plotfile', 'pdf', 'overview_24h.pdf'),
        arg('splice', 'cell line splice', arity='*',
            sprintf("splice/splice-%s_rbm24_vs_luc24.rds", c("H838", "H1650", "HCC70")))
    )

    diff_expr = readRDS(args$diff_expr) %>% mutate(cline = sub("_covar", "", cond), stype = "gex")
    splice_all = readRDS(args$splice_all) %>% mutate(cline = "all")
    splice_cline = sapply(args$splice, readRDS, simplify=FALSE)
    names(splice_cline) = sub(".*-([A-Z0-9]+)_.*", "\\1", names(splice_cline))
    splice_cline = bind_rows(splice_cline, .id="cline")

    splice = bind_rows(splice_all, splice_cline) %>%
        select(cline, everything()) %>%
        mutate(genes = lapply(genes, . %>% dplyr::rename(estimate=IncLevelDifference, statistic=stat, adj.p=FDR))) %>%
        tidyr::gather("collection", "res", -cline, -stype, -junction) %>%
        tidyr::unnest(res) %>%
        group_by(cline, stype, collection, label) %>%
            top_n(-1, adj.p) %>%
        ungroup()
    gex = diff_expr %>%
        filter(cond != "all") %>%
        mutate(junction = "covar",
               genes = lapply(genes, . %>% dplyr::rename(label=gene_name, estimate=log2FoldChange, statistic=stat, adj.p=padj))) %>%
        select(cline, stype, everything()) %>% select(-cond) %>%
        tidyr::gather("collection", "res", -cline, -stype, -junction) %>%
        tidyr::unnest(res)
    both = bind_rows(gex, splice) %>%
        mutate(cline = factor(cline, levels=c("all", "H838", "H1650", "HCC70")),
               stype = factor(stype, levels=c("gex", unique(splice$stype)))) %>%
        select(collection, cline, stype, junction, label, estimate, statistic, adj.p)

    use = both %>%
        filter(cline == "all", stype != "gex", adj.p < 0.05,
               # many similar complexes with similar names + separation stats (Nrp1 mouse?)
               !grepl("FARP2-NRP1-PlexinA[134] complex", label),
               !grepl("Nrp1", label),
               collection != "CORUM_all") %>%
        group_by(collection, stype) %>%
            top_n(-3, adj.p) %>%
        ungroup() %>%
        select(collection, label) %>%
        distinct()

    dset = inner_join(use, both) %>%
        mutate(statistic = sign(statistic) * pmin(abs(statistic), 5)) %>%
        group_by(cline, stype) %>%
            mutate(size = abs(statistic) / max(abs(statistic), na.rm=TRUE)) %>%
        ungroup()

    ggplot(dset, aes(x=stype, y=factor(label, levels=rev(unique(use$label))), fill=statistic)) +
        geom_point(aes(size=size, shape=junction), alpha=0.9, color="grey") +
        scale_shape_manual(values=c(covar=23, jc=21, jcec=22), guide=guide_legend(override.aes=list(size=3))) +
        scale_fill_distiller(palette="RdBu") +
        facet_grid(collection ~ cline, scales="free_y", space="free_y") +
        coord_cartesian(clip="off") +
        labs(size = "stat per\ncolumn", shape="type") +
        theme_classic() +
        theme(axis.text.x = element_text(angle=45, hjust=1),
              strip.text.y = element_text(angle=0),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
})
