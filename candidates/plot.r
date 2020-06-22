library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
util = import('./util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('d', 'dset', 'rds', '../merge/LUAD.rds'),
    opt('y', 'yaml', 'yaml', 'LUAD/top-genes.yaml'),
    opt('t', 'tissue', 'pan|TCGA identifier', 'LUAD'),
    opt('o', 'outfile', 'xlsx', '/dev/null'),
    opt('p', 'plotfile', 'pdf', 'LUAD/top-genes.pdf'))

select = yaml::read_yaml(args$yaml)
fits = select$methods
genes = select$genes #TODO: use right set if not only genes
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
    if (args$tissue == "pan") {
        cd$Name = NA # do not label cell lines in pan-can plots (too many)
        sizes = c(3, 1.5) # mut, wt
        alphas = c(0.5, 0.5)
    } else {
        cd = filter(cd, cohort == args$tissue) %>%
            mutate(Name = ifelse(expr == expr_orig & copies == copies_orig, Name,
                        sprintf("%s\n%.1f;%.1g", Name, copies_orig, expr_orig)))
        sizes = c(2, 3) # mut, wt
        alphas = c(0.8, 1)
    }
    abl = util$summary_ccle(cd, assocs, et=et, top=top)
    stats = assocs %>%
        filter(dset=="ccle", fit=="rlm3", cna=="amp") %>%
        transmute(gene=name, label=sprintf("%.0f%% comp R^2=%.2f", -estimate*100, rsq)) %>%
        inner_join(fracs %>% select(gene, min_reads)) %>%
        mutate(gene = factor(gene, levels=top))
    pccle = ggplot(cd, aes(x=copies, y=expr)) +
        annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, alpha=0.2, fill="yellow") +
        geom_vline(xintercept=2, color="grey") +
        geom_vline(xintercept=c(2-et,2+et,1+et,3-et), color="grey", linetype="dotted") +
        geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), size=1, linetype="dashed") +
        geom_point(aes(shape=mut, size=is.na(mut), alpha=is.na(mut), fill=meth_class), color="black") +
        ggrepel::geom_text_repel(aes(label=Name), size=1, alpha=0.5, segment.alpha=0.2) +
        geom_label(data=fracs, aes(x=x, y=max_reads, label=label, hjust=hjust),
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        geom_label(data=stats, aes(y=min_reads, label=label), x=2, hjust="center",
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        facet_wrap(~ gene, scales="free") +
        scale_color_manual(name="Compensation", guide="legend",
                           values=c("brown", "red", "blue"),
                           labels=c("full", "none", "observed (amp)")) +
        scale_fill_brewer(palette="RdBu", direction=-1,
                          labels=c("lowest", "low", "high", "highest")) +
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
    abl = util$summary_tcga(td, assocs, et=et, top=top)
    stats = assocs %>%
        filter(dset=="tcga", fit=="rlm3", cna=="amp", adj=="pur") %>%
        transmute(gene=name, label=sprintf("%.0f%% comp R^2=%.2f", -estimate*100, rsq)) %>%
        inner_join(fracs %>% select(gene, min_reads)) %>%
        mutate(gene = factor(gene, levels=top))
    ptcga = ggplot(td, aes(x=cancer_copies, y=expr)) +
        annotate("rect", xmin=2-et, xmax=2+et, ymin=-Inf, ymax=Inf, fill="#dedede") +
        stat_bin2d(aes(fill=after_stat(count)), bins=30) +
        geom_density2d(bins=20, color="#000000b0") +
        scale_fill_distiller(palette="Spectral", trans="log10") +
        geom_vline(xintercept=c(2-et,2+et,1+et,3-et), color="black", linetype="dotted") +
        geom_abline(data=abl, aes(intercept=intcp, slope=slope, color=type), size=1, linetype="dashed") +
        geom_point(data = td %>% filter(!is.na(mut)),
                   aes(shape=mut, size=is.na(mut)), color="black", alpha=1) +
        geom_label(data=fracs, aes(x=x, y=max_reads, label=label, hjust=hjust),
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        geom_label(data=stats, aes(y=min_reads, label=label), x=2, hjust="center",
                   label.size=NA, fill="#ffffff80", size=2.5, fontface="bold") +
        facet_wrap(~ gene, scales="free") +
        scale_color_manual(name="Compensation", guide="legend", na.value="#00000033",
                           values=c("brown", "red", "blue", "#000000ff"),
                           labels=c("full", "none", "observed (amp)", "x")) +
        scale_shape_manual(name="Mutation", guide="legend", na.value=21,
                           values=c(0, seq_along(levels(td$mut))[-1]),
                           labels=levels(td$mut)) +
        scale_size_manual(name="has mut", guide="none",
                           values=c(2, 1), labels=c("mut", "wt")) +
#        scale_alpha_continuous(trans="log", range=c(0.1, 0.5)) +
        labs(title = paste("cancer copy TCGA compensation;",
                           "98/99th% shown (expr/copies); dashed line model, solid capped data"),
             y = "normalized read count") +
        theme_classic()

    ###
    ### Methylation quantification
    ###
#    ct = bind_rows(ccle=cd, tcga=mutate(td, copies = cancer_copies), .id = "dset") %>%
#        select(cohort, gene, dset, purity, copies, expr, meth)
#    eu_cna = ct %>%
#        group_by(dset) %>%
#        mutate(facetx = "all", cna = case_when(
#            copies < 2-et ~ "del",
#            copies < 2+et ~ "eu",
#            copies >= 2+et ~ "amp",
#            TRUE ~ as.character(NA)), subs=cna) %>% ungroup() %>% na.omit()
#    cna_diff = eu_cna %>%
#        mutate(facetx = subs) %>%
#        filter(subs != "eu") %>%
#        group_by(dset, facetx) %>%
#        mutate(subs = case_when(
#            expr < median(expr, na.rm=TRUE) ~ "low", #TODO: fix for actual comp lines
#            expr > median(expr, na.rm=TRUE) ~ "high",
#            TRUE ~ as.character(NA))) %>% ungroup() %>% na.omit()
#    ct = bind_rows(eu_cna, cna_diff) %>%
#        mutate(facetx = factor(facetx, levels=c("del", "all", "amp")),
#               subs = factor(subs, levels=c("low", "high", "del", "eu", "amp")),
#               dset = factor(dset, levels=c("tcga", "ccle")))
#    pmeth = lapply(top, util$plot_meth_quant, ct=ct)

    ###
    ### save data underlying plots
    ###
#    writexl::write_xlsx(list(orf=orfdata, ccle=cd, tcga=td, meth=ct), args$outfile)
    writexl::write_xlsx(list(orf=orfdata, ccle=cd, tcga=td), args$outfile)

    ###
    ### actually plot
    ###
    print(patchwork::wrap_plots(lapply(overview, plt$try)) + plot_layout(guides="collect"))
    if (nrow(orfdata) > 0)
        print(porf)
    print(pccle)
    print(ptcga)
#    print(patchwork::wrap_plots(pmeth))
}

dev.off()
