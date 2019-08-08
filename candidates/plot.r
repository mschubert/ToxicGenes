library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('d', 'dset', 'rds', 'merge/pan/genes.rds'),
    opt('y', 'yaml', 'yaml', 'pan/top-genes.yaml'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('p', 'plotfile', 'pdf', 'pan/top-genes.pdf'))

quantile = function(x, ..., na.rm=TRUE) stats::quantile(x, ..., na.rm=na.rm)
select = yaml::read_yaml(args$yaml)
fits = select$methods
top = select$genes #TODO: use right set if not only genes
shapes = c("oe", "amp", "del", "all")
shape_i = c(21, 24, 25, 23)
et = yaml::read_yaml(args$config)$euploid_tol

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

#TODO: get the min of all top hits so we set limits?

dset = readRDS(args$dset) %>%
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
ccle_top = intersect(rownames(ccledata$copies), top)
cd = ccledata$clines %>%
    select(CCLE_ID, Name, Site_Primary, tcga_code) %>%
    left_join(reshape2::melt(ccledata$copies[ccle_top,], value.name="copies")) %>%
    left_join(reshape2::melt(ccledata$eset[ccle_top,], value.name="expr")) %>%
    mutate(expr = expr * copies/2, # undo normmatrix normalization
           gene = factor(gene, levels=top)) %>%
    group_by(gene) %>%
        filter(expr > quantile(expr, 0.05) & expr < quantile(expr, 0.95),
               copies > min(1, quantile(copies, 0.05)) & copies < max(3, quantile(copies, 0.95))) %>%
    ungroup()
if (args$tissue != "pan")
    cd = filter(cd, tcga_code == args$tissue)
abl = cd %>%
    group_by(gene) %>%
    summarize(med_expr = median(expr[abs(copies-2) < et]),
              none = 0.5 * med_expr,
              full = 0,
              observed = NA) %>%
    tidyr::gather("type", "slope", -gene, -med_expr) %>%
    mutate(intcp = ifelse(type == "full", med_expr, 0))
pccle = ggplot(cd, aes(x=copies, y=expr)) +
    annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, alpha=0.2, fill="yellow") +
    geom_vline(xintercept=2, color="grey") +
    geom_vline(xintercept=c(2-et,2+et), color="grey", linetype="dotted") +
    geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), linetype="dashed") +
    geom_point(alpha=0.3) +
    geom_smooth(aes(color="blue"), method="lm", color="blue") +
    facet_wrap(~ gene, scales="free") +
#    scale_fill_identity(name="CNA", guide="legend", labels="euploid") +
    scale_color_manual(name="Compensation", guide="legend",
                       values=c("brown", "red", "blue"),
                       labels=c("full", "none", "observed")) +
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
    transmute(sample = Var2, gene = Var1, meth = value) %>%
    group_by(gene) %>%
    mutate(meth = scale(meth)) %>%
#    mutate(meth = factor(cut_number(meth, n=5, labels=FALSE))) %>%
    ungroup() %>%
    mutate(meth = cut(meth, c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf)))
#levels(tcga_meth$meth) = c("lowest", "low", "normal", "high", "highest")
tcga_mut = tcga$mutations() %>%
    transmute(sample = Tumor_Sample_Barcode, gene = Hugo_Symbol, mut = factor(Variant_Classification))
names(dimnames(tcga_expr)) = c("gene", "sample")
names(dimnames(tcga_cns)) = c("gene", "sample")
td = reshape2::melt(tcga_expr, value.name="expr") %>%
    inner_join(reshape2::melt(tcga_cns, value.name="copies")) %>%
    inner_join(tcga$purity() %>% transmute(sample=Sample, purity=estimate)) %>%
    mutate(cancer_copies = (copies - 2) / purity + 2) %>%
    group_by(gene) %>%
        filter(expr > quantile(expr, 0.02) & expr < quantile(expr, 0.98)) %>%
        mutate(copies = ifelse(
                    copies > min(1, quantile(copies, 0.02)) & copies < max(3, quantile(copies, 0.98)),
                    copies, NA),
               cancer_copies = ifelse(
                    cancer_copies > min(1, quantile(cancer_copies, 0.02)) &
                    cancer_copies < max(3, quantile(cancer_copies, 0.98)),
                    cancer_copies, NA)) %>%
    ungroup() %>%
    left_join(tcga_mut) %>%
    left_join(tcga_meth)
abl = td %>%
    group_by(gene) %>%
    summarize(med_expr = median(expr[abs(cancer_copies-2) < et], na.rm=TRUE),
              none = 0.5 * med_expr,
              full = 0,
              observed = NA) %>%
    tidyr::gather("type", "slope", -gene, -med_expr) %>%
    mutate(intcp = ifelse(type == "full", med_expr, 0))
ptcga = ggplot(td, aes(x=cancer_copies, y=expr)) +
    annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, alpha=0.2, fill="yellow") +
    geom_vline(xintercept=2, color="grey") +
    geom_vline(xintercept=c(2-et,2+et), color="grey", linetype="dotted") +
    geom_point(data = td %>% filter(is.na(meth) | is.na(mut)),
               aes(shape=mut), fill="#ffffff00", color="#00000033") +
    geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), linetype="dashed") +
    geom_point(data = td %>% filter(!is.na(meth)), color="#ffffff00", shape=21,
               aes(fill=meth, alpha=meth)) +
    geom_point(data = td %>% filter(!is.na(mut)),
               aes(shape=mut, size=is.na(mut)), color="black", alpha=1) +
    geom_smooth(method="lm") +
    facet_wrap(~ gene, scales="free") +
    scale_color_manual(name="Compensation", guide="legend", na.value="#00000033",
                       values=c("brown", "red", "blue", "#000000ff"),
                       labels=c("full", "none", "observed", "x")) +
    scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                       values=c(0, seq_along(levels(td$mut))[-1]),
                       labels=levels(td$mut)) +
    scale_size_manual(name="has mut", guide="legend",
                       values=c(2, 1), labels=c("mut", "wt")) +
    scale_fill_brewer(name="z CpG meth", type="diverging", palette="RdBu",
                      direction=-1, na.value="#ffffff00") +
    scale_alpha_manual(guide="none", values = c(1,0.5,0,0.5,1)) +
    labs(title = paste("cancer copy TCGA compensation;",
                       "98th% shown (expr/copies); yellow=euploid"),
         y = "normalized read count")

###
### actually plot
###
pdf(args$plotfile, 16, 12) # patchworkGrob creates empty page in PDF
ov = overview # only way to get the legend to work
pg1 = patchworkGrob(
    ( ( ov[[1]] | ov[[2]] | ov[[3]] | ov[[4]] ) /
      ( ov[[5]] | ov[[6]] | ov[[7]] | ov[[8]] ) /
      ( ov[[9]] | ov[[10]] | ov[[11]] | ov[[12]] ) )
)
dev.off()
pdf(args$plotfile, 16, 12)
gridExtra::grid.arrange(pg1, ex_legend, ncol=2, widths=c(10,1))
print(porf)
print(pccle)
print(ptcga)
dev.off()
