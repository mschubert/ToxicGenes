library(dplyr)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')
idmap = import('process/idmap')

do_fit = function(gene, emat, copies, covar=NA) {
    on.exit(message("Error: ", gene))
    df = data.frame(expr = emat[gene,] / mean(emat[gene,], na.rm=TRUE) - 1,
                    erank = rank(emat[gene,]),
                    copies = copies[gene,],
                    covar = covar) %>%
        na.omit()

    has_covar = length(unique(na.omit(covar))) > 1
    if (has_covar) {
        f1 = erank ~ covar + copies
        f2 = expr ~ covar + copies
    } else {
        f1 = erank ~ copies
        f2 = expr ~ copies
    }

    prank = lm(f1, data=df) %>%
        broom::tidy() %>%
        filter(term == "copies") %>%
        pull(p.value)

    mod = MASS::rlm(f2, data=df, maxit=100) %>%
        broom::tidy() %>%
        filter(term == "copies") %>%
        select(-term) %>%
        mutate(n_aneup = sum(abs(df$copies-2) > 0.2),
               p.value = prank)

    on.exit()
    mod
}

all_fits = function(emat, copies, tissues=NA) {
    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
        dev = function(x) abs(x-2),
        all = identity
    )

    do_ffun = function(ffun) {
        copies2 = ffun(copies)
        fit_gene = function(g) do_fit(g, emat, copies2, tissues)
        tibble(gene = rownames(emat)) %>%
            mutate(res = parallel::mclapply(gene, fit_gene)) %>%
            tidyr::unnest() %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    }
    lapply(ffuns, do_ffun)
}

do_plot = function(data) {
    data %>%
        mutate(label=gene, size=n_aneup) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
        plt$volcano(base.size=0.2, label_top=50, repel=TRUE,
                    text.size=2.5, x_label_bias=5, pos_label_bias=0.15)
}

sys$run({
    args = sys$cmd$parse(
        opt('t', 'tissue', 'TCGA identifier', 'LUAD'),
        opt('o', 'outfile', 'xlsx', 'LUAD.xlsx'),
        opt('p', 'plotfile', 'pdf', 'LUAD.pdf'))

    # excl x,y chroms

    purity = tcga$purity() %>%
        filter(!is.na(estimate))

    reads = tcga$rna_seq(args$tissue) %>%
        tcga$filter(cancer=TRUE, primary=TRUE)
    rownames(reads) = idmap$gene(rownames(reads), to="hgnc_symbol")
    reads = reads[rowMeans(reads) >= 10,]

    if (args$tissue == "SKCM") #TODO:
        reads = reads[rownames(reads) != "TRGC1",]

    copies = tcga$cna_genes(args$tissue)
    rownames(copies) = idmap$gene(rownames(copies), to="hgnc_symbol")
    narray::intersect(purity$Sample, reads, copies, along=2)
#    copies[] = copies / narray::rrep(purity$estimate, nrow(copies)) #TODO: how to adjust this?
    copies[] = t(t(copies) / purity$estimate)
    copies = na.omit(copies)
    narray::intersect(reads, copies, along=1)

    cdata = data.frame(tissue = tcga$barcode2study(colnames(reads)))
    rownames(cdata) = colnames(reads)

    emat = DESeq2::DESeqDataSetFromMatrix(reads, cdata, ~1) %>%
        DESeq2::estimateSizeFactors(normMatrix=copies) %>%
        DESeq2::counts(normalized=TRUE)

#    if (args$tissue != "pan")
#        dset$idx$tcga_code[dset$idx$tcga_code != args$tissue] = NA

    fits = all_fits(emat, copies, purity$estimate)
    plots = lapply(fits, do_plot)

    pdf(args$plotfile)
    for (i in seq_along(plots))
        print(plots[[i]] + ggtitle(names(plots)[i]))
    dev.off()

    writexl::write_xlsx(fits, args$outfile)
})
