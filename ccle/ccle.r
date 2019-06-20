library(dplyr)
sys = import('sys')
plt = import('plot')

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
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('c', 'cores', 'parallel cores', getOption("mc.cores", 1L)),
        opt('t', 'tissue', 'TCGA identifier', 'pan'),
        opt('o', 'outfile', 'xlsx', 'pan.xlsx'),
        opt('p', 'plotfile', 'pdf', 'pan.pdf'))

    options(mc.cores = as.integer(args$cores))

    dset = readRDS(args$infile)
    emat = DESeq2::DESeqDataSetFromMatrix(dset$expr, dset$idx, ~1) %>%
        DESeq2::estimateSizeFactors(normMatrix=dset$copies) %>%
        DESeq2::counts(normalized=TRUE)

    if (args$tissue != "pan")
        dset$idx$tcga_code[dset$idx$tcga_code != args$tissue] = NA

    fits = all_fits(emat, dset$copies, dset$idx$tcga_code)
    plots = lapply(fits, do_plot)

    pdf(args$plotfile)
    for (i in seq_along(plots))
        print(plots[[i]] + ggtitle(names(plots)[i]))
    dev.off()

    writexl::write_xlsx(fits, args$outfile)
})
