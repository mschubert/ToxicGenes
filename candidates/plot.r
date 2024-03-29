library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
util = import('./util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('d', 'dset', 'rds', '../merge/pan.rds'),
    opt('y', 'yaml', 'yaml', 'manual.yaml'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('o', 'outfile', 'xlsx', '/dev/null'),
    opt('p', 'plotfile', 'pdf', 'pan/manual.pdf')
)

select = yaml::read_yaml(args$yaml)
fits = select$methods
genes = select$genes #TODO: use right set if not only genes
et = yaml::read_yaml(args$config)$euploid_tol
tcga_ccle_fit = setdiff(select$methods, "lm")
stopifnot(length(tcga_ccle_fit) == 1)
if (args$tissue == "COADREAD") {
    args$tissue = c("COAD", "READ")
} else if (args$tissue == "NSCLC") {
    args$tissue = c("LUAD", "LUSC")
}

comp_cols = c(full="brown", none="brown", amp="red", del="blue", all="darkviolet")

if (!is.list(genes))
    genes = list(genes=genes)

cairo_pdf(args$plotfile, 16, 12, onefile=TRUE)
for (i in seq_along(genes)) {
    top = genes[[i]]
    if (length(top) == 0)
        next
    print(plt$text(names(genes)[i], size=pmin(20, 100/sqrt(nchar(names(genes)[i])))))

    assocs = readRDS(args$dset)
    dset = assocs %>%
        filter(fit %in% fits) %>%
        mutate(statistic = pmax(statistic, -50),
               name = factor(name, levels=top))
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
        ggtitle("ORF drop-out (loess normalized, red line: mean +/- SD)") +
        theme_classic()

    ###
    ### CCLE data
    ###
    cd = util$load_ccle(top, et=et)
    fracs = util$frac_labels(cd, copies, et=et)
    ccle_assocs = assocs %>%
        filter(name %in% top, dset=="ccle", fit==tcga_ccle_fit) %>%
        mutate(name = factor(name, levels=top))
    if (args$tissue == "pan") {
        cd$Name = NA # do not label cell lines in pan-can plots (too many)
        sizes = c(3, 1.5) # mut, wt
        alphas = c(0.5, 0.5)
    } else {
        cd = filter(cd, cohort %in% args$tissue) %>%
            mutate(Name = ifelse(expr == expr_orig & copies == copies_orig, Name,
                        sprintf("%s\n%.1f;%.1g", Name, copies_orig, expr_orig)))
        sizes = c(2, 3) # mut, wt
        alphas = c(0.8, 1)
    }
    abl = util$comp_summary(ccle_assocs)
    stats = util$comp_stats(cd, copies, ccle_assocs, fracs)
    pccle = ggplot(cd, aes(x=copies, y=expr)) +
        annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, alpha=0.2, fill="yellow") +
        geom_vline(xintercept=2, color="grey") +
        geom_vline(xintercept=c(2-et,2+et,1+et,3-et), color="grey", linetype="dotted") +
        geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), size=1, linetype="dashed") +
        geom_point(aes(shape=mut, size=is.na(mut), alpha=is.na(mut), fill=meth), color="black") +
        ggrepel::geom_text_repel(aes(label=Name), size=1, alpha=0.5, segment.alpha=0.2) +
        geom_label(data=fracs, aes(x=x, y=max_reads, label=label, hjust=hjust),
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        geom_label(data=stats, aes(x=x, y=min_reads, label=label), hjust="center",
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        facet_wrap(~ gene, scales="free") +
        scale_color_manual(name="Compensation", guide="legend", values=comp_cols) +
        scale_fill_distiller(palette="RdBu", direction=-1) +
        guides(fill = guide_legend(override.aes=list(shape=21, size=5))) +
        scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                           values=c(0, seq_along(levels(cd$mut))[-1]),
                           labels=levels(cd$mut)) +
        scale_size_manual(guide="none", values=sizes) +
        scale_alpha_manual(guide="none", values=alphas) +
        labs(title = paste("CCLE compensation;",
                           "98/99th% shown (expr/copies); yellow=euploid"),
             y = "normalized read count") +
        theme_classic()

    ###
    ### TCGA data
    ###
    td = util$load_tcga(args$tissue, top, et=et)
    fracs = util$frac_labels(td, cancer_copies, et=et)
    tcga_assocs = assocs %>%
        filter(name %in% top, dset=="tcga", fit==tcga_ccle_fit, adj=="pur") %>%
        mutate(name = factor(name, levels=top))
    abl = util$comp_summary(tcga_assocs)
    stats = util$comp_stats(td, cancer_copies, tcga_assocs, fracs)
    ptcga = ggplot(td, aes(x=cancer_copies, y=expr, group=gene)) +
        annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, fill="#dedede") +
        stat_bin2d(aes(fill=after_stat(count)), bins=30) +
        geom_density2d(color="#000000b0", breaks=c(0.5,0.15,0.05), size=0.7, contour_var="ndensity") +
        scale_fill_distiller(palette="Spectral", trans="log10") +
        geom_vline(xintercept=c(2-et,2+et,1+et,3-et), color="black", linetype="dotted") +
        geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), size=1, linetype="dashed") +
        geom_point(data = td %>% filter(!is.na(mut)),
                   aes(shape=mut, size=is.na(mut)), color="black", alpha=1) +
        geom_label(data=fracs, aes(x=x, y=max_reads, label=label, hjust=hjust),
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        geom_label(data=stats, aes(x=x, y=min_reads, label=label), hjust="center",
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        facet_wrap(~ gene, scales="free") +
        scale_color_manual(name="Compensation", guide="legend", values=comp_cols) +
        scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                           values=c(0, seq_along(levels(td$mut))[-1]),
                           labels=levels(td$mut)) +
        scale_size_manual(name="has mut", guide="none",
                           values=c(2, 1), labels=c("mut", "wt")) +
#        scale_alpha_continuous(trans="log", range=c(0.1, 0.5)) +
        labs(title = "cancer copy TCGA compensation; 98/99th% shown (expr/copies)",
             y = "normalized read count") +
        theme_classic()

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
    pmeth = lapply(top, function(g) tryCatch(util$plot_meth_quant(ct, g),
                       error = function(e) plt$text(as.character(e))))

    ###
    ### save data underlying plots
    ###
    writexl::write_xlsx(list(orf=orfdata, ccle=cd, tcga=td, meth=ct), args$outfile)

    ###
    ### actually plot
    ###
    print(patchwork::wrap_plots(lapply(overview, plt$try)) + plot_layout(guides="collect"))
    if (nrow(orfdata) > 0)
        print(porf)
    print(pccle)
    print(ptcga)
    print(patchwork::wrap_plots(pmeth))
}

dev.off()
