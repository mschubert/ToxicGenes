library(dplyr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
sys = import('sys')
plt = import('plot')
util = import('./util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('d', 'dset', 'rds', 'merge/LUAD/genes.rds'),
    opt('y', 'yaml', 'yaml', 'LUAD/top-genes.yaml'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'LUAD'),
    opt('o', 'outfile', 'xlsx', '/dev/null'),
    opt('p', 'plotfile', 'pdf', 'LUAD/top-genes.pdf'))

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

    assocs = readRDS(args$dset)
    dset = assocs %>%
        filter(fit %in% fits) %>%
        mutate(statistic = pmax(statistic, -50),
               name = factor(name, levels=top))
    yaxis_floor = dset %>% filter(name %in% top) %>% pull(statistic) %>% min(na.rm=TRUE)
    overview = lapply(top, util$plot_stats, dset=dset)

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
    cd = util$load_ccle(top)
    if (args$tissue == "pan") {
        cd$Name = NA # do not label cell lines in pan-can plots (too many)
        sizes = c(3, 1.5) # mut, wt
        alphas = c(0.5, 0.5)
    } else {
        cd = filter(cd, cohort == args$tissue)
        sizes = c(2, 3) # mut, wt
        alphas = c(0.8, 1)
    }
    abl = util$summary_ccle(cd, assocs, et=et)
    pccle = ggplot(cd, aes(x=copies, y=expr)) +
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
    td = util$load_tcga(args$tissue, top, et=et)
    abl = util$summary_tcga(td, assocs, et=et)
    ptcga = ggplot(td, aes(x=cancer_copies, y=expr)) +
        util$stat_loess2d(aes(fill=meth_eup_scaled), se_size=TRUE) +
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
    pmeth = lapply(top, util$plot_meth_quant, ct=ct)

    ###
    ### save data underlying plots
    ###
    writexl::write_xlsx(list(orf=orfdata, ccle=cd, tcga=td, meth=ct), args$outfile)

    ###
    ### actually plot
    ###
    ex_legend = cowplot::get_legend(overview[[1]]) # only way to get the legend to work
    ov = lapply(seq_len(12+i-1), function(i) {
        if (i > length(overview))
            plot_spacer()
        else {
            p = overview[[i]]
            if (class(try(ggplot_build(p))) == "try-error")
                plot_spacer()
            else
                p + guides(color=FALSE, fill=FALSE, shape=FALSE)
        }
    })
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
