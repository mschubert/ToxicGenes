library(dplyr)
library(ggplot2)
library(patchwork)
idmap = import('process/idmap')
tcga = import('data/tcga')

quantile = function(x, ..., na.rm=TRUE) stats::quantile(x, ..., na.rm=na.rm)
shapes = c("amp"=24, "del"=25, "all"=23, "oe"=21)

#' Plot summary statistics of associations using FDR and percentiles
#'
#' @param dset  Association data set from merge/{tissue}/genes.rds
#' @param gene  Character string which gene to plot
#' @return      ggplot2 object
plot_stats = function(dset, gene) {
    yaxis_floor = dset %>% filter(name %in% gene) %>% pull(estimate) %>% min(na.rm=TRUE)
    cur = filter(dset, name == gene) %>%
        mutate(label = ifelse(adj %in% c("none", "pur"),
                              sprintf("%.2f th\nFDR %.1g", pctile, adj.p), NA))
    ggplot(dset, aes(x=1, y = estimate, color=adj)) +
        geom_hline(yintercept=0, linetype="dashed", color="grey") +
        geom_violin(position="identity", alpha=0) +
        geom_point(data=cur, aes(fill=adj, shape=cna),
                   color="black", size=3, alpha=0.5) +
        scale_shape_manual(values=shapes) +
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
load_ccle = function(top, et=0.15) {
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
                   expr = pmax(pmin(expr, quantile(expr, 0.98)), quantile(expr, 0.02)),
                   copies = pmin(copies, max(3+et, quantile(copies, 0.99))) %>%
                                 pmax(min(1-et, quantile(copies, 0.01)))) %>%
        ungroup()
}

#' Provide data for the compensation lines
#'
#' @param assocs  Associations for top genes in given data set
#' @return  Slopes and intercepts of lines
comp_summary = function(assocs) {
    dset = assocs %>%
        select(gene=name, type=cna, estimate, eup_reads)
    observed = dset %>%
        mutate(slope = 0.5 * eup_reads * (1 + estimate),
               intcp = eup_reads - 2*slope)
    ref = dset %>%
        group_by(gene) %>%
        summarize(eup_reads = median(eup_reads, na.rm=TRUE)) %>%
        mutate(none = 0.5 * eup_reads,
               full = 0) %>%
        tidyr::gather("type", "slope", -gene, -eup_reads) %>%
        mutate(intcp = eup_reads - 2*slope)
    bind_rows(observed, ref)
}

#' Construct a label with compensation statistics
#'
#' @param dset  data.frame with underlying data
#' @param dset_x  field in 'dset' that is plotted on x axis
#' @param assocs  Associations for top genes in given data set
#' @param fracs  data.frame containing min_reads for positioning
#' @return  data.frame for label and its positions
comp_stats = function(dset, dset_x, assocs, fracs) {
    dx = dset %>%
        group_by(gene) %>%
        summarize(x = mean(range({{ dset_x }}, na.rm=TRUE)))
    if (!"rsq" %in% colnames(assocs)) # currently not present in stan-nb
        assocs$rsq = NA
    assocs %>%
        select(name, cna, estimate, rsq) %>%
        mutate(estimate = -100*estimate) %>%
        tidyr::pivot_wider(names_from=cna, values_from=c(estimate, rsq)) %>%
#        transmute(gene=name, label=sprintf("comp/R^2 ▲ %.0f%% %.2f ▼ %.0f%% %.2f ◆ %.0f%% %.2f",
#               estimate_amp, rsq_amp, estimate_del, rsq_del, estimate_all, rsq_all)) %>%
        transmute(gene=name, label=sprintf("comp ▲ %.0f%%",
               estimate_amp, rsq_amp)) %>%
        inner_join(fracs %>% select(gene, min_reads)) %>%
        inner_join(dx)
}

#' Label fractions of amplified/euploid/deleted samples
#'
#' @param dset    A data.frame for candidate plotting
#' @param copies  Field of x axis
frac_labels = function(dset, copies, et=0.15) {
    fracs = dset %>%
        group_by(gene) %>%
        mutate(ntot = n(),
               min_reads = min(expr, na.rm=TRUE),
               max_reads = max(expr, na.rm=TRUE),
               cuts = cut({{ copies }}, c(-Inf, 1+et, 2-et, 2+et, 3-et, Inf), labels=FALSE),
               x = setNames(c(1+et, 1.5, 2, 2.5, 3-et), 1:5)[as.character(cuts)],
               hjust = setNames(c("right", rep("center", 3), "left"), 1:5)[as.character(cuts)]) %>%
        group_by(gene, x, hjust, min_reads, max_reads) %>%
            summarize(n=n(), pct=n()/ntot[1]) %>%
        ungroup() %>%
        mutate(label = sprintf("%i\n%.0f%%", n, pct*100))
}

#' Load TCGA data
#'
#' @param cohort  TCGA cohort identifier or "pan"
#' @param genes   Character vector of genes of interest
load_tcga = function(cohort, top, et=0.15) {
    if (cohort == "pan") cohort2 = tcga$cohorts() else cohort2 = cohort
    load_expr = function(cohort, genes) {
        mat = lapply(cohort2, tcga$rna_seq) %>%
            narray::stack(along=2, allow_overwrite=TRUE) %>%
            DESeq2::DESeqDataSetFromMatrix(data.frame(id=colnames(.)), ~1) %>%
                DESeq2::estimateSizeFactors() %>%
                DESeq2::counts(normalized=TRUE)
        mat = mat[rowSums(mat) != 0,] # gene needs to have at least 1 read
        rownames(mat) = idmap$gene(rownames(mat), to="hgnc_symbol")
        mat[intersect(genes, rownames(mat)),,drop=FALSE]
    }
    load_copies = function(cohort, genes) {
        lapply(cohort2, function(x) {
            cns = tcga$cna_genes(x, gene="external_gene_name")
            cns[intersect(genes, rownames(cns)),,drop=FALSE]
        }) %>% narray::stack(along=2)
    }
    mut_proc = function(cohort) {
        tcga$mutations(cohort) %>%
            transmute(sample = substr(Tumor_Sample_Barcode, 1, 16),
                      gene = Hugo_Symbol,
                      mut = factor(Variant_Classification))
    }
    tcga_expr = load_expr(cohort, top) # only primary because purity() only
    tcga_cns = load_copies(cohort, top) # has those samples
    tcga_meth = tcga$meth(cohort2, cpg="avg", mvalues=TRUE, genes=top) %>% # cpg=stdev too many NAs
        reshape2::melt() %>%
        transmute(sample = Var2, gene = Var1, meth = value)
    tcga_mut = lapply(cohort2, mut_proc) %>% bind_rows()
    names(dimnames(tcga_expr)) = c("gene", "sample")
    names(dimnames(tcga_cns)) = c("gene", "sample")

    surv = tcga$clinical() %>%
        transmute(sample = paste0(submitter_id, "-01A"),
                  os_status = factor(vital_status, levels=c("alive", "dead")),
                  os_days = pmax(days_to_last_follow_up, days_to_death, na.rm=TRUE))

    scale_ref = function(x, ref) {
        medref = median(x[ref], na.rm=TRUE)
        sdref = sd(x[ref], na.rm=TRUE)
        (x - medref) / sdref
    }

    subtype = tcga$immune() %>%
        transmute(sample=paste0(barcode, "-01A"), subtype=`TCGA Subtype`) %>%
        na.omit()

    td = reshape2::melt(tcga_expr, value.name="expr") %>%
        mutate(cohort = tcga$barcode2study(sample)) %>%
        inner_join(reshape2::melt(tcga_cns, value.name="copies")) %>%
        inner_join(tcga$purity() %>% transmute(sample=Sample, purity=consensus)) %>%
        mutate(expr = expr / purity,
               cancer_copies = (copies - 2) / purity + 2) %>%
        group_by(gene) %>%
#            filter(expr > quantile(expr, 0.01) & expr < quantile(expr, 0.99)) %>%
            mutate(expr = pmax(pmin(expr, quantile(expr, 0.98)), quantile(expr, 0.02)),
                   copies = pmin(copies, max(3+et, quantile(copies, 0.99))) %>%
                                 pmax(min(1-et, quantile(copies, 0.01))),
                   cancer_copies = pmin(cancer_copies, max(3+et, quantile(cancer_copies, 0.99))) %>%
                                 pmax(min(1-et, quantile(cancer_copies, 0.01)))) %>%
        ungroup() %>%
        left_join(tcga_mut) %>%
        left_join(tcga_mut %>% filter(gene == "TP53") %>% transmute(sample=sample, p53_mut=mut)) %>%
        left_join(tcga_meth) %>%
        left_join(surv) %>%
        left_join(subtype) %>%
        group_by(gene, cohort) %>%
            mutate(meth_eup_scaled = scale_ref(meth, abs(cancer_copies-2)<et),
                   meth_eup_scaled = pmax(pmin(meth_eup_scaled, 2), -2)) %>%
        ungroup() %>%
        mutate(gene = factor(gene, levels=top))
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
#' @param gamma     Parameter to control smoothness of the surface (higher: more smooth)
#' @param se_alpha  Alpha of tiles corresponds to confidence (default: FALSE)
#' @param se_size   Size of tiles corresponds to confidence (default: FALSE)
#' @param cap_z     Only allow gradient as extreme as data (default: TRUE)
#' @param by        Variable to condition on (prediction will be on 1)
#' @param ...       Passed to ggplot layer
stat_gam2d = function(mapping = NULL, data = NULL, geom = "tile",
        position = "identity", na.rm = FALSE, show.legend = NA, gamma=10,
        inherit.aes = TRUE, bins = 20, se_alpha=FALSE, se_size=FALSE, cap_z=TRUE, ...) {
    gam2d = ggproto("gam2d", Stat, # also: sd/se as alpha?; pts<=density?
        compute_group = function(data, scales, gamma=10, bins=20, se_alpha=FALSE, se_size=FALSE, cap_z=TRUE) {
            rx = range(data$x, na.rm=TRUE)
            ry = range(data$y, na.rm=TRUE)
            df = expand.grid(x = seq(rx[1], rx[2], length.out=bins),
                             y = seq(ry[1], ry[2], length.out=bins))

            if ("by" %in% colnames(data)) {
                lsurf = mgcv::gam(fill ~ s(x, y, by=by), data=data, gamma=gamma)
                df$by = 1
            } else {
                lsurf = mgcv::gam(fill ~ s(x, y), data=data, gamma=gamma)
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
        params = list(bins = bins, gamma=gamma, se_alpha=se_alpha, se_size=se_size, cap_z=cap_z,
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
