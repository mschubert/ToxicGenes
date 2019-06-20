library(dplyr)
io = import('io')
sys = import('sys')
plt = import('plot')

do_fit = function(gene) {
    on.exit(message("Error: ", gene))
    df = data.frame(expr = nmat[gene,] / mean(nmat[gene,], na.rm=TRUE) - 1,
                    erank = rank(nmat[gene,]),
                    copies = dset$copies[gene,],
                    tcga = dset$idx$tcga_code) %>%
        na.omit()

    prank = lm(erank ~ tcga + copies, data=df) %>%
        broom::tidy() %>%
        filter(term == "copies") %>%
        pull(p.value)

    mod = MASS::rlm(expr ~ tcga + copies, data=df, maxit=100) %>%
        broom::tidy() %>%
        filter(term == "copies") %>%
        select(-term) %>%
        mutate(size = sum(df$copies > 2.2),
               p.value = prank)

    on.exit()
    mod
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', '../data/ccle/dset.rds'),
        opt('o', 'outfile', 'xlsx', 'pan.xlsx'),
        opt('p', 'plotfile', 'pdf', 'pan.pdf'))

    # amp/del/abs/all?

    dset = readRDS(args$infile)
    emat = DESeq2::DESeqDataSetFromMatrix(dset$expr, dset$idx, ~1) %>%
        DESeq2::estimateSizeFactors(normMatrix=dset$copies) %>%
        DESeq2::counts(normalized=TRUE)
    nmat = emat
    nmat[dset$copies < 1.8] = NA # amps only

    pancov = tibble(gene = rownames(nmat)) %>%
        mutate(res = purrr::map(gene, do_fit)) %>%
        tidyr::unnest() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)

    p = pancov %>%
        mutate(label=gene) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
        plt$volcano(base.size=0.2, label_top=50, repel=TRUE,
                    x_label_bias=5, pos_label_bias=0.15)

    pdf(args$plotfile)
    print(p)
    dev.off()

    writexl::write_xlsx(pancov, args$outfile)
})
