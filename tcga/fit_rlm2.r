library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')
gset = import('data/genesets')

models = function(type, covar) {
    fmls = list(
        naive = list(
            "TRUE" = expr - cancer_copies ~ covar + cancer_copies,
            "FALSE" = expr - cancer_copies ~ cancer_copies
        ),
        pur = list (
            "TRUE" = expr - cancer_copies ~ covar + cancer + cancer_copies,
            "FALSE" = expr - cancer_copies ~ cancer + cancer_copies
        ),
        puradj = list( # 0 + cancer + stroma is singular?!
            "TRUE" = expr - cancer_copies ~ covar + stroma + cancer_copies,
            "FALSE" = expr - cancer_copies ~ stroma + cancer_copies
        )
    )
    fmls[[type]][[as.character(covar)]]
}

do_fit = function(genes, fml, emat, copies, purity, covar=0, et=0.15) {
    df = data.frame(gene = rep(genes, each=ncol(emat)),
                    expr = c(emat[genes,,drop=FALSE]),
                    cancer_copies = c((copies[genes,] - 2) / purity + 2),
                    purity = rep(purity, length(genes)),
                    covar = rep(covar, length(genes))) %>%
        na.omit() %>%
        mutate(cancer = purity,
               stroma = 1 - purity)

    tryCatch({
        df2 = df %>%
            filter(cancer_copies > 2-et & cancer_copies < 2+et) %>%
            group_by(gene, covar) %>%
            filter(n() >= 5) %>% # at least 5 euploid samples per cohort per gene
            summarize(intcp = MASS::rlm(expr ~ stroma, maxit=100)$coefficients["(Intercept)"]) %>%
            ungroup() %>%
            inner_join(df, by=c("gene", "covar")) %>%
            mutate(expr = expr / intcp - 1,
                   cancer_copies = cancer_copies / 2 - 1) %>%
            na.omit() %>%
            sample_n(min(nrow(.), 1e5))

        mobj = MASS::rlm(fml, data=df2, maxit=100)
        mod = broom::tidy(mobj) %>%
            filter(term == "cancer_copies") %>%
            select(-term) %>%
            mutate(n_aneup = sum(abs(df$cancer_copies-2) > et),
                   n_genes = length(genes),
                   p.value = sfsmisc::f.robftest(mobj, var="cancer_copies")$p.value)
    }, error = function(e) {
        warning(conditionMessage(e), immediate.=TRUE)
        data.frame(estimate = NA)
    })
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('t', 'tissue', 'TCGA identifier', 'LUAD'),
        opt('y', 'type', 'naive|pur|puradj', 'naive'),
        opt('s', 'setfile', 'rds', '../data/genesets/CH.HALLMARK.rds'),
        opt('j', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '4096'),
        opt('o', 'outfile', 'xlsx', 'LUAD/genes.xlsx'))

    et = yaml::read_yaml(args$config)$euploid_tol

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

    if (grepl("genes\\.xlsx", args$outfile))
        sets = setNames(rownames(emat), rownames(emat))
    else
        sets = readRDS(args$setfile) %>%
            gset$filter(min=4, valid=rownames(emat))

    w = clustermq::workers(n_jobs = min(as.integer(args$cores),
                                        ceiling(length(sets)/20)),
                           template = list(memory = as.integer(args$memory)))

    ffuns = list(
        amp = function(x) { x[x < 2-et] = NA; x },
        del = function(x) { x[x > 2+et] = NA; x },
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        has_covar = length(unique(cdata$tissue)) != 1
        res = clustermq::Q(do_fit, genes=sets, workers=w, pkgs="dplyr",
                const = list(fml=models(args$type, has_covar),
                             emat=emat, copies=ff(copies), et=et,
                             purity=purity$estimate, covar=cdata$tissue)) %>%
            setNames(names(sets)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
    w$cleanup()
})
