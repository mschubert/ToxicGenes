library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')
gset = import('data/genesets')

models = function(type, covar) {
    fmls = list(
        naive = list(
            "TRUE" = erank ~ 0 + covar + cancer_copy_dev,
            "FALSE" = erank ~ 0 + cancer_copy_dev
        ),
        pur = list (
            "TRUE" = erank ~ 0 + covar + cancer + cancer_copy_dev,
            "FALSE" = erank ~ 0 + cancer + cancer_copy_dev
        ),
        puradj = list(
            "TRUE" = erank ~ 0 + covar + stroma + cancer + cancer_copy_dev,
            "FALSE" = erank ~ 0 + stroma + cancer + cancer_copy_dev
        )
    )
    fmls[[type]][[as.character(covar)]]
}

do_fit = function(genes, fml, emat, copies, purity, covar=0) {
    df = data.frame(expr = c(emat[genes,]),
                    cancer_copies = c((copies[genes,] - 2) / purity + 2),
                    purity = rep(purity, length(genes)),
                    covar = rep(covar, length(genes))) %>%
        na.omit() %>%
        sample_n(min(nrow(.), 1e5)) %>%
        mutate(expr = expr / cancer_copies,
               stroma = 2 * (1 - purity) / cancer_copies,
               cancer = (purity * cancer_copies) / cancer_copies, # simplifies to CCF
               cancer_copy_dev = cancer_copies - 2,
               erank = rank(expr) / nrow(.) - 0.5)

    mod = lm(fml, data=df) %>%
        broom::tidy() %>%
        filter(term == "cancer_copy_dev") %>%
        select(-term) %>%
        mutate(n_aneup = sum(abs(df$cancer_copies-2) > 0.2),
               n_genes = length(genes))
}

sys$run({
    args = sys$cmd$parse(
        opt('t', 'tissue', 'TCGA identifier', 'LUAD'),
        opt('y', 'type', 'naive|pur|puradj', 'naive'),
        opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
        opt('c', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '4096'),
        opt('o', 'outfile', 'xlsx', 'LUAD/genes.xlsx'))

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
#    emat = emat / rowMeans(emat, na.rm=TRUE) - 1

    emat = narray::map(emat, along=2, subsets=cdata$tissue,
                       function(x) rank(x) / length(x) - 0.5)

    if (grepl("genes\\.xlsx", args$outfile))
        sets = setNames(rownames(emat), rownames(emat))
    else
        sets = readRDS(args$setfile) %>%
            gset$filter(min=4, valid=rownames(emat))

    w = clustermq::workers(n_jobs = min(as.integer(args$cores),
                                        ceiling(length(sets)/20)),
                           template = list(memory = as.integer(args$memory)))

    ffuns = list(
        amp = function(x) { x[x < 1.8] = NA; x },
        del = function(x) { x[x > 2.2] = NA; x },
#        dev = function(x) abs(x-2),
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        has_covar = length(unique(cdata$tissue)) != 1
        res = clustermq::Q(do_fit, genes=sets, workers=w, pkgs="dplyr",
                const = list(fml=models(args$type, has_covar),
                             emat=emat, copies=ff(copies),
                             purity=purity$estimate, covar=cdata$tissue)) %>%
            setNames(names(sets)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
    w$cleanup()
})
