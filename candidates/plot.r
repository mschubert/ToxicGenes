library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('d', 'dset', 'rds', 'merge/pan/genes.rds'),
    opt('y', 'yaml', 'yaml', 'pan/top-genes.yaml'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('p', 'plotfile', 'pdf', 'pan/top-genes.pdf'))

quantile = function(x, ..., na.rm=TRUE) stats::quantile(x, ..., na.rm=na.rm)
select = yaml::read_yaml(args$yaml)
fits = select$methods
genes = select$genes #TODO: use right set if not only genes
shapes = c("oe", "amp", "del", "all")
shape_i = c(21, 24, 25, 23)
et = yaml::read_yaml(args$config)$euploid_tol

if (!is.list(genes))
    genes = list(genes=genes)

pdf(args$plotfile, 16, 12)
for (i in seq_along(genes)) {
top = genes[[i]]
print(plt$text(names(genes)[i], size=20))

#' Get the percentile of x in y
plot_stats = function(gene) {
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
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
}

assocs = readRDS(args$dset)
dset = assocs %>%
    filter(fit %in% fits) %>%
    mutate(statistic = pmax(statistic, -50),
           name = factor(name, levels=top))
yaxis_floor = dset %>% filter(name %in% top) %>% pull(statistic) %>% min(na.rm=TRUE)
overview = lapply(top, plot_stats)
ex_legend = cowplot::get_legend(overview[[1]])
overview = lapply(overview, function(p) p + guides(color=FALSE, fill=FALSE, shape=FALSE))

###
### ORF data
###
orfdata = readRDS("../orf/overview.rds")
if (args$tissue != "pan")
    orfdata = filter(orfdata, tissue == args$tissue)
orfdata = orfdata %>%
    filter(`GENE SYMBOL` %in% top) %>%
    mutate(`GENE SYMBOL` = factor(`GENE SYMBOL`, levels=top)) %>%
    group_by(`GENE SYMBOL`) %>%
        mutate(construct_i = as.integer(factor(`Construct IDs`))) %>%
    ungroup()
porf = ggplot(orfdata, aes(x=DMSO, y=z_LFC)) +
    geom_hline(yintercept=0, color="red") +
    geom_hline(yintercept=c(-1,1), color="red", linetype="dotted") +
    geom_point(aes(color=tissue, shape=factor(construct_i)), size=3, alpha=0.6) +
    facet_wrap(~ `GENE SYMBOL`, drop=FALSE) +
    ggtitle("ORF drop-out (loess normalized, red line: mean +/- SD)")

###
### CCLE data
###
ccledata = readRDS("../data/ccle/dset.rds")
names(dimnames(ccledata$copies)) = c("gene", "CCLE_ID")
names(dimnames(ccledata$eset)) = c("gene", "CCLE_ID")
names(dimnames(ccledata$meth)) = c("gene", "CCLE_ID")
cd = ccledata$clines %>%
    select(CCLE_ID, Name, Site_Primary, cohort=tcga_code) %>%
    left_join(ccledata$copies[intersect(top, rownames(ccledata$copies)),] %>%
              reshape2::melt(value.name="copies")) %>%
    left_join(ccledata$eset[intersect(top, rownames(ccledata$eset)),] %>%
              reshape2::melt(value.name="expr")) %>%
    left_join(ccledata$meth[intersect(top, rownames(ccledata$meth)),] %>%
              reshape2::melt(value.name="meth")) %>%
    left_join(ccledata$mut %>% dplyr::rename(Name=cline, mut=type)) %>%
    mutate(expr = expr * copies/2, # undo normmatrix normalization
           gene = factor(gene, levels=top),
           mut = factor(mut),
           purity = 1) %>%
    group_by(gene) %>%
        mutate(expr = pmax(pmin(expr, quantile(expr, 0.95)), quantile(expr, 0.05)),
               copies = pmax(pmin(copies, quantile(copies, 0.95)), quantile(copies, 0.05))) %>%
    ungroup() %>%
    group_by(gene) %>%
        mutate(meth_class = rank(meth, ties="min", na="keep") / sum(!is.na(meth)),
               meth_class = cut(meth_class, c(0, 0.25, 0.5, 0.75, 1), labels=FALSE),
               meth_class = factor(meth_class)) %>%
    ungroup()
if (args$tissue == "pan") {
    cd$Name = NA # do not label cell lines in pan-can plots (too many)
    sizes = c(3, 1.5) # mut, wt
    alphas = c(0.5, 0.5)
} else {
    cd = filter(cd, cohort == args$tissue)
    sizes = c(2, 3) # mut, wt
    alphas = c(0.8, 1)
}
abl = cd %>%
    group_by(gene) %>%
    summarize(med_expr = median(expr[abs(copies-2) < et]),
              none = 0.5 * med_expr,
              full = 0) %>%
    left_join(assocs %>% filter(dset=="ccle", fit=="rlm2", cna=="all") %>%
              transmute(gene=factor(name, levels=top), observed=estimate)) %>%
    mutate(observed = none * (1 + observed)) %>%
    tidyr::gather("type", "slope", -gene, -med_expr) %>%
    mutate(intcp = med_expr - 2*slope)
pccle =
    ggplot(cd, aes(x=copies, y=expr)) +
    annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, alpha=0.2, fill="yellow") +
    geom_vline(xintercept=2, color="grey") +
    geom_vline(xintercept=c(2-et,2+et), color="grey", linetype="dotted") +
    geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), size=1, linetype="dashed") +
    geom_point(aes(shape=mut, size=is.na(mut), alpha=is.na(mut), fill=meth_class), color="black") +
    ggrepel::geom_text_repel(aes(label=Name), size=1, alpha=0.5, segment.alpha=0.2) +
    facet_wrap(~ gene, scales="free") +
    scale_color_manual(name="Compensation", guide="legend",
                       values=c("brown", "red", "blue"),
                       labels=c("full", "none", "observed")) +
    scale_fill_brewer(palette="RdBu", direction=-1,
                      labels=c("lowest", "low", "high", "highest")) +
    guides(fill = guide_legend(override.aes=list(shape=21, size=5))) +
    scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                       values=c(0, seq_along(levels(cd$mut))[-1]),
                       labels=levels(cd$mut)) +
    scale_size_manual(guide="none", values=sizes) +
    scale_alpha_manual(guide="none", values=alphas) +
    labs(title = paste("CCLE compensation;",
                       "95th% shown (expr/copies); yellow=euploid"),
         y = "normalized read count")

###
### TCGA data
###
idmap = import('process/idmap')
tcga = import('data/tcga')
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
tcga_expr = load_expr(args$tissue, top) # only primary because purity() only
tcga_cns = load_copies(args$tissue, top) # has those samples
tcga_meth = tcga$meth(tcga$cohorts(), cpg="avg", mvalues=TRUE, genes=top) %>% # cpg=stdev too many NAs
    reshape2::melt() %>%
    transmute(cohort = tcga$barcode2study(Var2), sample = Var2, gene = Var1, meth = value) %>%
    group_by(gene) %>%
        mutate(meth_class = scale(meth)) %>%
#    group_by(cohort, gene) %>%
#    mutate(meth = factor(cut_number(meth, n=5, labels=FALSE))) %>%
    ungroup() %>%
    mutate(meth_class = cut(meth, c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf)))
#levels(tcga_meth$meth) = c("lowest", "low", "normal", "high", "highest")
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
    left_join(tcga_meth) %>%
    group_by(gene) %>%
        mutate(meth_eup_scaled = scale_ref(meth, abs(cancer_copies-2)<et),
               meth_eup_scaled = pmax(pmin(meth_eup_scaled, 2), -2)) %>%
    ungroup() %>%
    mutate(gene = factor(gene, levels=top))
abl = td %>%
    group_by(gene) %>%
    summarize(med_expr = median(expr[abs(cancer_copies-2) < et], na.rm=TRUE),
              none = 0.5 * med_expr,
              full = 0) %>%
    left_join(assocs %>% filter(dset=="tcga", fit=="rlm2", cna=="all") %>%
              transmute(gene=factor(name, levels=top), observed=estimate)) %>%
    mutate(observed = none * (1 + observed)) %>%
    tidyr::gather("type", "slope", -gene, -med_expr) %>%
    mutate(intcp = med_expr - 2*slope)
stat_loess2d = function(mapping = NULL, data = NULL, geom = "tile",
        position = "identity", na.rm = FALSE, show.legend = NA,
        inherit.aes = TRUE, bins = 20, cap_z=TRUE, ...) {
    loess2d = ggproto("loess2d", Stat, # also: sd/se as alpha?; pts<=density?
        compute_group = function(data, scales, bins=20, cap_z=TRUE) {
            rx = range(data$x, na.rm=TRUE)
            ry = range(data$y, na.rm=TRUE)
            df = expand.grid(x = seq(rx[1], rx[2], length.out=bins),
                             y = seq(ry[1], ry[2], length.out=bins))

            lsurf = loess(fill ~ x + y, data=data)
            df$fill = c(predict(lsurf, newdata=df))
            if (cap_z) {
                df$fill = pmax(df$fill, min(df$fill, na.rm=TRUE))
                df$fill = pmin(df$fill, max(df$fill, na.rm=TRUE))
            }

            df$width = rep(diff(rx) / (bins - 1), bins)
            df$height = rep(diff(ry) / (bins - 1), each=bins)
            df
        },
        required_aes = c("x", "y", "fill")
    )

    layer(
        stat = loess2d, data = data, mapping = mapping, geom = geom,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(bins = bins, cap_z=TRUE, na.rm = na.rm, ...)
    )
}
ptcga =
    ggplot(td, aes(x=cancer_copies, y=expr)) +
    stat_loess2d(aes(fill=meth_eup_scaled)) +
    geom_density2d(bins=20, color="#00000050") +
    geom_vline(xintercept=c(2-et,2+et), color="black", linetype="dotted") +
    geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), size=1, linetype="dashed") +
    geom_point(data = td %>% filter(!is.na(mut)),
               aes(shape=mut, size=is.na(mut)), color="black", alpha=1) +
    facet_wrap(~ gene, scales="free") +
    scale_fill_gradient2(low="blue", mid="white", high="red") +
    scale_color_manual(name="Compensation", guide="legend", na.value="#00000033",
                       values=c("brown", "red", "blue", "#000000ff"),
                       labels=c("full", "none", "observed", "x")) +
    scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                       values=c(0, seq_along(levels(td$mut))[-1]),
                       labels=levels(td$mut)) +
    scale_size_manual(name="has mut", guide="none",
                       values=c(2, 1), labels=c("mut", "wt")) +
    scale_alpha_continuous(trans="log", range=c(0.1, 0.5)) +
    labs(title = paste("cancer copy TCGA compensation;",
                       "98th% shown (expr/copies); dashed line model, solid capped data"),
         y = "normalized read count")

###
### Methylation quantification
###
plot_gene = function(g) {
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
        ggtitle(g)
}
ct = bind_rows(ccle=cd, tcga=mutate(td, copies = cancer_copies), .id = "dset") %>%
    select(cohort, gene, dset, purity, copies, expr, meth)
eu_cna = ct %>%
    group_by(dset) %>%
    mutate(facetx = "all", cna = case_when(
        copies < 2-et ~ "del",
        copies < 2+et ~ "eu",
        copies >= 2+et ~ "amp",
        TRUE ~ as.character(NA)), subs=cna) %>% ungroup() %>% na.omit()
cna_diff = eu_cna %>%
    mutate(facetx = subs) %>%
    filter(subs != "eu") %>%
    group_by(dset, facetx) %>%
    mutate(subs = case_when(
        expr < median(expr, na.rm=TRUE) ~ "low", #TODO: fix for actual comp lines
        expr > median(expr, na.rm=TRUE) ~ "high",
        TRUE ~ as.character(NA))) %>% ungroup() %>% na.omit()
ct = bind_rows(eu_cna, cna_diff) %>%
    mutate(facetx = factor(facetx, levels=c("del", "all", "amp")),
           subs = factor(subs, levels=c("low", "high", "del", "eu", "amp")),
           dset = factor(dset, levels=c("tcga", "ccle")))
pmeth = lapply(top, plot_gene)

###
### actually plot
###
ov = overview # only way to get the legend to work
pdf("/dev/null")
pg1 = patchworkGrob(
    ( ( ov[[1]] | ov[[2]] | ov[[3]] | ov[[4]] ) /
      ( ov[[5]] | ov[[6]] | ov[[7]] | ov[[8]] ) /
      ( ov[[9]] | ov[[10]] | ov[[11]] | ov[[12]] ) )
)
dev.off()
gridExtra::grid.arrange(pg1, ex_legend, ncol=2, widths=c(10,1))
print(porf)
print(pccle)
print(ptcga)
print(cowplot::plot_grid(plotlist=pmeth))

}

dev.off()
