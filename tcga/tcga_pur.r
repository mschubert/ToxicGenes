library(dplyr)
library(patchwork)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')
idmap = import('process/idmap')

do_fit = function(gene, emat, copies, purity, covar=0) {
    on.exit(message("Error: ", gene))
    df = data.frame(expr = emat[gene,] / mean(emat[gene,], na.rm=TRUE),
                    purity = purity,
                    cancer_copies = (copies[gene,] - 2) / purity + 2,
                    covar = covar) %>%
        na.omit() %>%
        mutate(expr = expr / cancer_copies,
               stroma = 2 * (1 - purity) / cancer_copies,
               cancer = (purity * cancer_copies) / cancer_copies) # simplifies to CCF

    n_aneup = sum(abs(df$cancer_copies-2) > 0.2)
    if (n_aneup < 5) {
        on.exit(message("Not enough aneuploid samples: ", gene))
        return(data.frame(estimate=NA))
    }

    if (length(unique(na.omit(covar))) > 1) {
        fml = expr ~ covar + cancer + cancer_copies #+ delta_cancer + delta_stroma
    } else {
        fml = expr ~ cancer + cancer_copies #+ stroma:cancer_copies
    }
    mobj = MASS::rlm(fml, data=df, maxit=100)
    mod = broom::tidy(mobj) %>%
        mutate(p.value = NA)
    mod$p.value[mod$term == "cancer_copies"] = sfsmisc::f.robftest(mobj, var="cancer_copies")$p.value
#    mod$p.value[mod$term == "stroma"] = sfsmisc::f.robftest(mobj, var="stroma")$p.value

    mod = mod %>%
        filter(term == "cancer_copies") %>%
        select(-term) %>%
        mutate(n_aneup = n_aneup)

    on.exit()
    mod
}

all_fits = function(emat, copies, purity, tissues=0) {
    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
#        dev = function(x) abs(x-2),
        all = identity
    )

    do_ffun = function(ffun) {
        copies2 = ffun(copies)
        fit_gene = function(g) do_fit(g, emat, copies2, purity$estimate, tissues)
        tibble(gene = rownames(emat)) %>%
            mutate(res = lapply(gene, fit_gene)) %>%
            tidyr::unnest() %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    }
    lapply(ffuns, do_ffun)
}

do_plot = function(data, title) {
    df = data %>%
        mutate(label=gene, size=n_aneup) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate")
    p1 = df %>%
        filter(estimate < 0) %>%
        plt$volcano(base.size=0.2, label_top=30, repel=TRUE) +
        labs(title = sprintf("%s (%i aneup genes)", title, sum(!is.na(df$estimate))),
             subtitle = "compensated")
    p2 = df %>%
        filter(estimate > 0) %>%
        plt$volcano(base.size=0.2, label_top=30, repel=TRUE) +
        labs(subtitle="hyperactivated") +
        theme(axis.title.y = element_blank())
    p1 + p2 + plot_layout(nrow=1)
}

sys$run({
    args = sys$cmd$parse(
        opt('t', 'tissue', 'TCGA identifier', 'LUAD'),
        opt('o', 'outfile', 'xlsx', 'LUAD.xlsx'),
        opt('p', 'plotfile', 'pdf', 'LUAD.pdf'))

    if (args$tissue == "pan")
        args$tissue = tcga$cohorts()

    # excl x,y chroms

    purity = tcga$purity() %>%
        filter(!is.na(estimate))

    reads = lapply(args$tissue, tcga$rna_seq) %>%
        narray::stack(along=2) %>%
        tcga$filter(cancer=TRUE, primary=TRUE)
    rownames(reads) = idmap$gene(rownames(reads), to="hgnc_symbol")
    reads = reads[rowMeans(reads) >= 10,]

    copies = lapply(args$tissue, tcga$cna_genes) %>%
        narray::stack(along=2)
    rownames(copies) = idmap$gene(rownames(copies), to="hgnc_symbol")
    narray::intersect(purity$Sample, reads, copies, along=2)
#    copies[] = copies / narray::rrep(purity$estimate, nrow(copies)) #TODO: how to adjust this?
#    copies[] = t(t(copies) / purity$estimate)
    copies = na.omit(copies)
    narray::intersect(reads, copies, along=1)

    cdata = data.frame(tissue = tcga$barcode2study(colnames(reads)))
    rownames(cdata) = colnames(reads)

    emat = DESeq2::DESeqDataSetFromMatrix(reads, cdata, ~1) %>%
        DESeq2::estimateSizeFactors() %>% # total ploidy to scale lib size
        DESeq2::counts(normalized=TRUE)

    fits = all_fits(emat, copies, purity)
#    fits2 = unlist(lapply(fits, function(x) split(x, x$term)), recursive=FALSE)

    pdf(args$plotfile, 10, 8)
    for (i in seq_along(fits))
        print(do_plot(fits[[i]], names(fits)[i]))
    dev.off()

    writexl::write_xlsx(fits, args$outfile)
})
