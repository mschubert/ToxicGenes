import_package("dplyr", attach=TRUE)
import_package("ggplot2", attach=TRUE)
util = import('../candidates/util')

#' Handle error as warning and return NULL
muffle = function(e) { warning(e, immediate.=TRUE); NULL }

#' Get number of cases in matrix
nc = function(mat) { re = ncol(mat); if (is.null(re)) 0 else re }

plot_l2d = function(dset, variable, et=0.15, from=NA, to=NA, by="purity") {
    # https://slowkow.com/notes/ggplot2-color-by-density/
    get_density = function(x, y, ...) {
        dens = MASS::kde2d(x, y, ...)
        ix = findInterval(x, dens$x)
        iy = findInterval(y, dens$y)
        ii = cbind(ix, iy)
        dens$z[ii]
    }
    lowdens = dset %>%
        select(sample, cohort, mut, p53_mut, cancer_copies, expr) %>%
        filter(!is.na(cancer_copies) & !is.na(expr)) %>%
        group_by(cohort, p53_mut) %>%
            mutate(dens = get_density(cancer_copies, expr)) %>%
            filter(dens < quantile(dens, 0.25) | !is.na(mut)) %>%
        ungroup()

    if (all(na.omit(dset[[variable]]) >= 0)) {
        fill = scale_fill_viridis_c(option="magma", direction=-1, limits=c(from, to))
    } else {
        fill = scale_fill_gradientn(colours=rev(RColorBrewer::brewer.pal(11,"RdBu")),
                                limits=c(from, to))
    }
    ggplot(dset, aes(x=cancer_copies, y=expr)) +
        util$stat_gam2d(aes_string(fill=variable, by=by), se_size=TRUE) +
        geom_density2d(breaks=c(0.5,0.15,0.05), color="chartreuse4", size=0.7, contour_var="ndensity") +
        geom_vline(xintercept=c(2-et,2+et), color="springgreen4", linetype="dotted", size=1.5) +
        facet_grid(p53_mut ~ cohort, scales="free") +
        fill +
        geom_point(data=lowdens, aes(shape=mut), color="magenta", alpha=0.6, size=3) +
        scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                           values=c(0, seq_along(levels(dset$mut))[-1]),
                           labels=levels(dset$mut)) +
        guides(fill=guide_legend(title="")) +
        labs(title = variable,
             y = "normalized read count") +
        theme_classic()
}
