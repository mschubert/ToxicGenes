library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

models = function(type, covar) {
    fmls = list(
        naive = list(
            "TRUE" = expr ~ covar + cancer_copies,
            "FALSE" = expr ~ cancer_copies
        ),
        pur = list (
            "TRUE" = expr ~ covar + cancer + cancer_copies,
            "FALSE" = expr ~ cancer + cancer_copies
        ),
        puradj = list(
            "TRUE" = expr ~ covar + stroma + cancer + cancer_copies,
            "FALSE" = expr ~ stroma + cancer + cancer_copies
        )
    )
    fmls[[type]][[as.character(covar)]]
}

fit_gene = function(gene, fml, emat, copies, purity, covar=0) {
    on.exit(message("Error: ", gene))
    df = data.frame(expr = emat[gene,], purity = purity, covar = covar,
                    cancer_copies = (copies[gene,] - 2) / purity + 2) %>%
        na.omit() %>%
        mutate(expr = expr / cancer_copies,
               stroma = 2 * (1 - purity) / cancer_copies,
               cancer = (purity * cancer_copies) / cancer_copies) # simplifies to CCF

    mobj = MASS::rlm(fml, data=df, maxit=100)
    mod = broom::tidy(mobj) %>%
        filter(term == "cancer_copies") %>%
        select(-term) %>%
        mutate(n_aneup = sum(abs(df$cancer_copies-2) > 0.2),
               size = n_aneup,
               p.value = sfsmisc::f.robftest(mobj, var="cancer_copies")$p.value)

    on.exit()
    mod
}

sys$run({
    args = sys$cmd$parse(
        opt('t', 'tissue', 'TCGA identifier', 'LUAD'),
        opt('y', 'type', 'naive|pur|puradj', 'naive'),
        opt('c', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '4096'),
        opt('o', 'outfile', 'xlsx', 'LUAD.xlsx'))

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
    copies = na.omit(copies[!is.na(rownames(copies)),])
    narray::intersect(reads, copies, along=1)

    cdata = data.frame(tissue = tcga$barcode2study(colnames(reads)))
    rownames(cdata) = colnames(reads)

    emat = DESeq2::DESeqDataSetFromMatrix(reads, cdata, ~1) %>%
        DESeq2::estimateSizeFactors() %>% # total ploidy to scale lib size
        DESeq2::counts(normalized=TRUE)
    emat = emat / rowMeans(emat, na.rm=TRUE) - 1

    w = clustermq::workers(n_jobs = as.integer(args$cores),
                           template = list(memory = as.integer(args$memory)))

    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
#        dev = function(x) abs(x-2),
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        cmat = ff(copies)
        fml = models(args$type, length(unique(cdata$tissue)) != 1)

        res = tibble(gene = rownames(emat)) %>%
            mutate(res = clustermq::Q(fit_gene, gene=gene,
                workers=w, pkgs="dplyr",
                const = list(fml=fml, emat=emat, copies=cmat,
                             purity=purity$estimate, covar=cdata$tissue))) %>%
            tidyr::unnest() %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
})
