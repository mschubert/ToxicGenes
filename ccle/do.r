library(dplyr)
io = import('io')
sys = import('sys')
plt = import('plot')

dset = readRDS("../data/ccle/dset.rds")

emat = DESeq2::DESeqDataSetFromMatrix(dset$expr, dset$idx, ~1) %>%
    DESeq2::estimateSizeFactors(normMatrix=dset$copies) %>%
    DESeq2::counts(normalized=TRUE)

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

nmat = emat
nmat[dset$copies < 1.8] = NA # amps only

pancov = tibble(gene = rownames(nmat)) %>%
    mutate(res = purrr::map(gene, do_fit)) %>%
    tidyr::unnest() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

x = pancov %>%
    mutate(label=gene) %>%
    plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
    plt$volcano(base.size=0.2, label_top=50, repel=TRUE,
                x_label_bias=5, pos_label_bias=0.15)

#pdf(args$outfile)
print(x)
dev.off()

# writexl::write_ ...
