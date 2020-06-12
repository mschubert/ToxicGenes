library(dplyr)
library(ggplot2)
library(patchwork)
idmap = import('process/idmap')
tcga = import('data/tcga')

quantile = function(x, ..., na.rm=TRUE) stats::quantile(x, ..., na.rm=na.rm)

#' Plot summary statistics of associations using FDR and percentiles
#'
#' @param dset  Association data set from merge/{tissue}/genes.rds
#' @param gene  Character string which gene to plot
#' @return      ggplot2 object
plot_stats = function(dset, gene) {
    cur = filter(dset, name == gene) %>%
        mutate(cna = factor(cna, levels=shapes),
               label = ifelse(adj %in% c("none", "puradj"),
                              sprintf("%.2f th\nFDR %.1g", pctile, adj.p), NA))
    ggplot(dset, aes(x=1, y = statistic, color=adj)) +
        geom_hline(yintercept=0, linetype="dashed", color="grey") +
        geom_violin(position="identity", alpha=0) +
        geom_point(data=cur, aes(fill=adj, shape=cna),
                   color="black", size=3, alpha=0.5) +
        scale_shape_manual(values=shape_i) +
        ggrepel::geom_text_repel(data=cur, size=2, box.padding=unit(7, "pt"),
            aes(label=label), color="black", direction="y", segment.alpha=0.3) +
        facet_wrap(~ dset + fit, scale="free_x", nrow=1) +
        coord_cartesian(ylim=c(0,yaxis_floor)) +
        labs(title = gene) +
        theme_classic() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
}

#' Load CCLE data
#'
#' @param top  Character vector of genes to load
load_ccle = function(top) {
    ccledata = readRDS("../data/ccle/dset.rds")
    names(dimnames(ccledata$copies)) = c("gene", "CCLE_ID")
    names(dimnames(ccledata$eset)) = c("gene", "CCLE_ID")
    names(dimnames(ccledata$meth)) = c("gene", "CCLE_ID")
    p53mut = ccledata$mut %>%
        filter(gene == "TP53") %>%
        transmute(Name=cline, p53_mut=type) %>%
        na.omit()
    cd = ccledata$clines %>%
        select(CCLE_ID, Name, Site_Primary, cohort=tcga_code) %>%
        left_join(ccledata$copies[intersect(top, rownames(ccledata$copies)),] %>%
                  reshape2::melt(value.name="copies")) %>%
        left_join(ccledata$eset[intersect(top, rownames(ccledata$eset)),] %>%
                  reshape2::melt(value.name="expr")) %>%
        left_join(ccledata$meth[intersect(top, rownames(ccledata$meth)),] %>%
                  reshape2::melt(value.name="meth")) %>%
        left_join(ccledata$mut %>% dplyr::rename(Name=cline, mut=type)) %>%
        left_join(p53mut) %>%
        mutate(expr = expr * copies/2, # undo normmatrix normalization
               gene = factor(gene, levels=top),
               mut = factor(mut),
               purity = 1) %>%
        group_by(gene) %>%
            mutate(expr_orig = expr, copies_orig = copies,
                   expr = pmax(pmin(expr, quantile(expr, 0.95)), quantile(expr, 0.05)),
                   copies = pmax(pmin(copies, quantile(copies, 0.95)), quantile(copies, 0.05))) %>%
        ungroup() %>%
        group_by(gene) %>%
            mutate(meth_class = rank(meth, ties="min", na="keep") / sum(!is.na(meth)),
                   meth_class = cut(meth_class, c(0, 0.25, 0.5, 0.75, 1), labels=FALSE),
                   meth_class = factor(meth_class)) %>%
        ungroup()
}

summary_ccle = function(cd, assocs, et=0.15) {
    to_merge = assocs %>%
        filter(dset=="ccle", fit=="rlm3", cna=="amp") %>%
        transmute(gene=factor(name, levels=top), observed=estimate, eup_reads=eup_reads)
    abl = cd %>%
        select(gene) %>%
        left_join(to_merge) %>%
        mutate(none = 0.5 * eup_reads,
               full = 0,
               observed = none * (1 + observed)) %>%
        tidyr::gather("type", "slope", -gene, -eup_reads) %>%
        mutate(intcp = eup_reads - 2*slope)
}

#' Load TCGA data
#'
#' @param cohort  TCGA cohort identifier or "pan"
#' @param genes   Character vector of genes of interest
load_tcga = function(cohort, top, et=0.15) {
    load_expr = function(cohort, genes) {
        if (cohort == "pan")
            cohort = tcga$cohorts()
        lapply(cohort, function(x) {
            expr = tcga$rna_seq(x)
            rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
            expr[intersect(genes, rownames(expr)),,drop=FALSE]
        }) %>% narray::stack(along=2)
    }
    load_copies = function(cohort, genes) {
        if (cohort == "pan")
            cohort = tcga$cohorts()
        lapply(cohort, function(x) {
            cns = tcga$cna_genes(x, gene="external_gene_name")
            cns[intersect(genes, rownames(cns)),,drop=FALSE]
        }) %>% narray::stack(along=2)
    }
    tcga_expr = load_expr(cohort, top) # only primary because purity() only
    tcga_cns = load_copies(cohort, top) # has those samples
    if (cohort == "pan")
        tcga_meth = tcga$meth(tcga$cohorts(), cpg="avg", mvalues=TRUE, genes=top) # cpg=stdev too many NAs
    else
        tcga_meth = tcga$meth(cohort, cpg="avg", mvalues=TRUE, genes=top) # cpg=stdev too many NAs
    tcga_meth = tcga_meth %>%
        reshape2::melt() %>%
        transmute(sample = Var2, gene = Var1, meth = value)
    tcga_mut = tcga$mutations() %>%
        transmute(sample = Tumor_Sample_Barcode, gene = Hugo_Symbol, mut = factor(Variant_Classification))
    names(dimnames(tcga_expr)) = c("gene", "sample")
    names(dimnames(tcga_cns)) = c("gene", "sample")

    scale_ref = function(x, ref) {
        medref = median(x[ref], na.rm=TRUE)
        sdref = sd(x[ref], na.rm=TRUE)
        (x - medref) / sdref
    }

    td = reshape2::melt(tcga_expr, value.name="expr") %>%
        mutate(cohort = tcga$barcode2study(sample)) %>%
        inner_join(reshape2::melt(tcga_cns, value.name="copies")) %>%
        inner_join(tcga$purity() %>% transmute(sample=Sample, purity=estimate)) %>%
        mutate(cancer_copies = (copies - 2) / purity + 2) %>%
        group_by(gene) %>%
            filter(expr > quantile(expr, 0.02) & expr < quantile(expr, 0.98)) %>%
            mutate(expr = pmax(pmin(expr, quantile(expr, 0.98)), quantile(expr, 0.02)),
                   copies = pmax(pmin(copies, quantile(copies, 0.98)), quantile(copies, 0.02)),
                   cancer_copies = pmax(pmin(copies, quantile(cancer_copies, 0.98)),
                                        quantile(cancer_copies, 0.02))) %>%
        ungroup() %>%
        left_join(tcga_mut) %>%
        left_join(tcga_mut %>% filter(gene == "TP53") %>% transmute(sample=sample, p53_mut=mut)) %>%
        left_join(tcga_meth) %>%
        group_by(gene, cohort) %>%
            mutate(meth_eup_scaled = scale_ref(meth, abs(cancer_copies-2)<et),
                   meth_eup_scaled = pmax(pmin(meth_eup_scaled, 2), -2)) %>%
        ungroup() %>%
        mutate(gene = factor(gene, levels=top))
}

#' TCGA expected vs. observed expression
#'
#' @param td  data.frame from load_tcga()
#' @return    data.frame for lines of expected vs observed compensation
summary_tcga = function(td, assocs, et=0.15) {
    to_merge = assocs %>%
        filter(dset=="tcga", fit=="rlm3", cna=="amp", adj=="puradj") %>%
        transmute(gene=factor(name, levels=top), observed=estimate, eup_reads=eup_reads)
    abl = td %>%
        select(gene) %>%
        left_join(to_merge) %>%
        mutate(none = 0.5 * eup_reads,
               full = 0,
               observed = none * (1 + observed)) %>%
        tidyr::gather("type", "slope", -gene, -eup_reads) %>%
        mutate(intcp = eup_reads - 2*slope)
}

#' 2D loess plot showing gradient on a third variable
#'
#' @param bins      Number of bins to partition data in x and y axis (default: 20)
#' @param se_alpha  Alpha of tiles corresponds to confidence (default: FALSE)
#' @param se_size   Size of tiles corresponds to confidence (default: FALSE)
#' @param cap_z     Only allow gradient as extreme as data (default: TRUE)
stat_loess2d = function(mapping = NULL, data = NULL, geom = "tile",
        position = "identity", na.rm = FALSE, show.legend = NA,
        inherit.aes = TRUE, bins = 20, se_alpha=FALSE, se_size=FALSE, cap_z=TRUE, ...) {
    loess2d = ggproto("loess2d", Stat, # also: sd/se as alpha?; pts<=density?
        compute_group = function(data, scales, bins=20, se_alpha=FALSE, se_size=FALSE, cap_z=TRUE) {
            rx = range(data$x, na.rm=TRUE)
            ry = range(data$y, na.rm=TRUE)
            df = expand.grid(x = seq(rx[1], rx[2], length.out=bins),
                             y = seq(ry[1], ry[2], length.out=bins))

            lsurf = loess(fill ~ x + y, data=data, span=0.1)
            pred = predict(lsurf, newdata=df, se=TRUE)
            df$fill = c(pred$fit)
            cor_factor = c(pred$se.fit) / 1.5
            cor_factor = pmax(0.1, 1 - abs(cor_factor/df$fill))

            df$width = rep(diff(rx) / (bins - 1), bins)
            df$height = rep(diff(ry) / (bins - 1), each=bins)

            if (cap_z) {
                df$fill = pmax(df$fill, min(data$fill, na.rm=TRUE))
                df$fill = pmin(df$fill, max(data$fill, na.rm=TRUE))
            }
            if (se_alpha)
                df$alpha = cor_factor
            if (se_size) {
                df$width = df$width * cor_factor
                df$height = df$height * cor_factor
            }
            df
        },
        required_aes = c("x", "y", "fill")
    )

    layer(
        stat = loess2d, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(bins = bins, se_alpha=se_alpha, se_size=se_size, cap_z=cap_z,
                      na.rm = na.rm, ...)
    )
}

#' 2D GAM plot showing gradient on a third variable
#'
#' @param bins      Number of bins to partition data in x and y axis (default: 20)
#' @param se_alpha  Alpha of tiles corresponds to confidence (default: FALSE)
#' @param se_size   Size of tiles corresponds to confidence (default: FALSE)
#' @param cap_z     Only allow gradient as extreme as data (default: TRUE)
#' @param by        Variable to condition on (prediction will be on 1)
stat_gam2d = function(mapping = NULL, data = NULL, geom = "tile",
        position = "identity", na.rm = FALSE, show.legend = NA,
        inherit.aes = TRUE, bins = 20, se_alpha=FALSE, se_size=FALSE, cap_z=TRUE, ...) {
    gam2d = ggproto("gam2d", Stat, # also: sd/se as alpha?; pts<=density?
        compute_group = function(data, scales, bins=20, se_alpha=FALSE, se_size=FALSE, cap_z=TRUE) {
            rx = range(data$x, na.rm=TRUE)
            ry = range(data$y, na.rm=TRUE)
            df = expand.grid(x = seq(rx[1], rx[2], length.out=bins),
                             y = seq(ry[1], ry[2], length.out=bins))

            if ("by" %in% colnames(data)) {
                lsurf = mgcv::gam(fill ~ s(x, y, by=by), data=data, select=TRUE)
                df$by = 1
            } else {
                lsurf = mgcv::gam(fill ~ s(x, y), data=data, select=TRUE)
            }
            pred = predict(lsurf, newdata=df, se=TRUE)
            df$fill = c(pred$fit)
            eff_max = quantile(data$fill, 0.98, na.rm=TRUE)
            eff_min = quantile(data$fill, 0.02, na.rm=TRUE)
            cor_factor = pmax(0.1, 1 - 2*(c(pred$se.fit)/(eff_max-eff_min)))

            df$width = rep(diff(rx) / (bins - 1), bins)
            df$height = rep(diff(ry) / (bins - 1), each=bins)

            if (cap_z) {
                df$fill = pmin(df$fill, eff_max)
                df$fill = pmax(df$fill, eff_min)
            }
            if (se_alpha)
                df$alpha = cor_factor
            if (se_size) {
                df$width = df$width * cor_factor
                df$height = df$height * cor_factor
            }
            df
        },
        required_aes = c("x", "y", "fill", "by") #TODO: by as optional aes
    )

    layer(
        stat = gam2d, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(bins = bins, se_alpha=se_alpha, se_size=se_size, cap_z=cap_z,
                      na.rm = na.rm, ...)
    )
}

#' 2D tile convolution for survival fit
#'
#' @param bins      Number of bins to partition data in x and y axis (default: 20)
#' @param se_alpha  Alpha of tiles corresponds to confidence (default: FALSE)
#' @param se_size   Size of tiles corresponds to confidence (default: FALSE)
stat_surv2d = function(mapping = NULL, data = NULL, geom = "tile",
        position = "identity", na.rm = FALSE, show.legend = NA,
        inherit.aes = TRUE, bins = 20, se_alpha=FALSE, se_size=FALSE, ...) {
    surv2d = ggproto("surv2d", Stat, # also: sd/se as alpha?; pts<=density?
        compute_group = function(data, scales, bins=20, se_alpha=FALSE, se_size=FALSE) {
            rx = range(data$x, na.rm=TRUE)
            ry = range(data$y, na.rm=TRUE)
            df = expand.grid(x = seq(rx[1], rx[2], length.out=bins),
                             y = seq(ry[1], ry[2], length.out=bins))
            df$width = rep(diff(rx) / (bins - 1), bins)
            df$height = rep(diff(ry) / (bins - 1), each=bins)

            compute_gridpt = function(i, j) {
#                keep_x = with(data, x > )
#                keep_y = with(data, y ...)
                res = survival::coxph(survival::Surv(time, status) ~ 1,
                                      data=data[keep_x & keep_y])
                broom::tidy(res)
            }

            lsurf = mgcv::gam(fill ~ s(x, y), data=data, select=TRUE)
            pred = predict(lsurf, newdata=df, se=TRUE)
            df$fill = c(pred$fit)

            cor_factor = pmax(0.1, 1 - 2*(c(pred$se.fit)/(eff_max-eff_min)))

            if (se_alpha)
                df$alpha = cor_factor
            if (se_size) {
                df$width = df$width * cor_factor
                df$height = df$height * cor_factor
            }
            df
        },
        required_aes = c("x", "y", "status", "time")
    )

    layer(
        stat = surv2d, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(bins = bins, se_alpha=se_alpha, se_size=se_size, cap_z=cap_z,
                      na.rm = na.rm, ...)
    )
}

#' Quantify methylation differences
#'
#' @param ct  Merged CCLE and TCGA data.frame
#' @return    ggplot2 object
plot_meth_quant = function(ct, g) {
    #TODO: lines pancan + per tissue
    cur = ct %>% filter(gene == g) %>%
        group_by(dset, gene) %>%
        mutate(meth = pmax(meth, quantile(meth, 0.1)),
               meth = pmin(meth, quantile(meth, 0.9)))
    stats = tibble(facetx = factor(c("del", "amp", "all", "all"), levels=levels(cur$facetx)),
                   xref = factor(c("low", "low", "eu", "eu"), levels=levels(cur$subs)),
                   xcmp = factor(c("high", "high", "del", "amp"), levels=levels(cur$subs))) %>%
        inner_join(cur, by="facetx") %>%
        filter(subs == xref | subs == xcmp) %>%
        group_by(dset, facetx, xref, xcmp) %>%
        summarize(res = list(tryCatch(error = function(e) tibble(intcp=NA, slope=NA, p.value=NA), {
            if (length(unique(cohort)) == 1)
                mod = broom::tidy(lm(meth ~ (1-purity) + (subs==xcmp)))
            else
                mod = broom::tidy(lm(meth ~ (1-purity) + cohort + (subs==xcmp)))
            tibble(intcp = mod$estimate[mod$term == "(Intercept)"],
                   slope = mod$estimate[mod$term == "subs == xcmpTRUE"],
                   p.value = mod$p.value[mod$term == "subs == xcmpTRUE"]) })
        )) %>%
        tidyr::unnest() %>%
        inner_join(cur %>% group_by(dset) %>% summarize(yrange = diff(range(meth, na.rm=TRUE)))) %>%
        mutate(label_col = ifelse(!is.na(p.value) & p.value < 0.05, "p<0.05", "ns"),
               label = sprintf("%.2g", p.value),
               angle = 0)
#               angle = atan(yrange/4 * slope / (xcmp-xref) * (180/pi)))
    ggplot(cur, aes(x=subs, y=meth), color="grey50") +
#        ggbeeswarm::geom_quasirandom() + # too many points
        geom_violin(aes(fill=cna)) +
        geom_boxplot(width=0.25, outlier.shape=NA) +
        geom_segment(data=stats, aes(x=xref, y=intcp, xend=xcmp, yend=intcp+slope,
                                     color=label_col), size=1) +
        geom_text(data=stats, size=2, aes(label=label, color=label_col, angle=angle,
                  # https://github.com/tidyverse/ggplot2/issues/577
                  x=ifelse(facetx == "all", (as.integer(xref)+as.integer(xcmp))/2 - 2,
                                            (as.integer(xref)+as.integer(xcmp))/2),
                  y=intcp+slope/2), vjust=-0.5) +
        facet_grid(dset ~ facetx, scales="free", space="free_x") +
        scale_color_manual(labels=c("ns", "p<0.05"), values=c("grey20", "blue"), guide=FALSE) +
        scale_fill_manual(labels=c("amp", "del", "eu"), guide=FALSE,
                          values=c("#83242455", "#3A3A9855", "white")) +
        ggtitle(g) +
        theme_classic()
}
