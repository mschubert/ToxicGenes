library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')
gset = import('data/genesets')

do_fit = function(genes, emat, copies, purity, covar=0, et=0.15, type="pur") {
    df2 = data.frame(gene = rep(genes, each=ncol(emat)),
                    expr = c(emat[genes,,drop=FALSE]),
                    cancer_copies = c((copies[genes,] - 2) / purity + 2),
                    purity = rep(purity, length(genes)),
                    covar = rep(covar, length(genes))) %>%
        na.omit() %>%
        mutate(covar = paste(covar, gene, sep=":"),
               cancer = purity,
               stroma = 1 - purity)

    tryCatch({
        reps = replicate(100, simplify=FALSE, {
            df = df2[sample(seq_len(nrow(df2)), replace=TRUE),]
            if (length(unique(na.omit(covar))) > 1) {
                fml = switch(type,
                    naive = expr ~ covar + cancer_copies,
                    pur = expr ~ covar + purity + cancer_copies,
                    puradj = expr ~ covar + purity + cancer_copies
                )
                fml_mean = switch(type,
                    naive = expr ~ covar,
                    pur = expr ~ covar + purity,
                    puradj = expr ~ covar + purity
                )
                mobj = lm(fml, data=df)
                ucovar = unique(df$covar)
                # should I predict purity=1 here or whatever purity there is in a sample?
                # does it matter? becausse we're regressing out purity again anyway
                pred = predict(mobj, newdata=expand.grid(covar=ucovar, purity=1, cancer_copies=2))
                expr_per_copy = data.frame(covar=ucovar, expr_per_copy=pred/2)
                df = inner_join(df, expr_per_copy, by="covar")
            } else {
                fml = switch(type,
                    naive = expr ~ cancer_copies,
                    pur = expr ~ purity + cancer_copies,
                    puradj = expr ~ purity + cancer_copies
                )
                fml_mean = switch(type,
                    naive = expr ~ 1,
                    pur = expr ~ purity,
                    puradj = expr ~ purity
                )
                mobj = lm(fml, data=df)
                # should I predict purity=1 here or whatever purity there is in a sample?
                # does it matter? becausse we're regressing out purity again anyway
                pred = predict(mobj, newdata=data.frame(purity=1, cancer_copies=2))
                df$expr_per_copy = pred / 2
            }

            if (mean(pred) < 0) # maybe add check for mobj2/intercept as well?
                stop("predicted negative expression for euploid")

            df$expr = with(df, expr - cancer_copies * expr_per_copy)
            mobj2 = lm(fml, data=df)
            mmean = lm(fml_mean, data=df)
            res = broom::tidy(mobj2) %>%
                filter(term == "cancer_copies") %>%
                select(-term) %>%
                mutate(estimate = 2 * estimate / mean(pred), # pct_comp
                       n_aneup = sum(abs(df$cancer_copies-2) > et),
                       n_genes = length(genes),
                       eup_reads = mean(pred),
                       slope_diff = 2 * estimate,
                       rsq = broom::glance(mobj2)$r.squared,
                       p.value = p.value)
        })

        res = bind_rows(reps) %>% arrange(statistic)
        res = res[round(nrow(res)/2),]
    }, error = function(e) {
        warning(genes, ": ", conditionMessage(e), immediate.=TRUE)
        data.frame(estimate = NA)
    })
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('t', 'tissue', 'TCGA identifier', 'LUAD'),
        opt('y', 'type', 'naive|pur|puradj', 'naive'),
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

    genes = setNames(rownames(emat), rownames(emat))

    w = clustermq::workers(n_jobs = 100, log=T, #_file="cmq.%a.log", #as.integer(args$cores),
                           template = list(memory = 40960)) #as.integer(args$memory)))

    ffuns = list(
        amp = function(x) { x[x < 2-et] = NA; x },
        del = function(x) { x[x > 2+et] = NA; x },
        all = identity
    )
    fits = lapply(ffuns, function(ff) {
        res = clustermq::Q(do_fit, genes=genes, workers=w, pkgs="dplyr",
                const = list(emat=emat, copies=ff(copies), et=et,
                             purity=purity$estimate, covar=cdata$tissue,
                             type=args$type)) %>%
            setNames(names(genes)) %>%
            bind_rows(.id="name") %>%
            mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
            arrange(adj.p, p.value)
    })

    writexl::write_xlsx(fits, args$outfile)
    w$cleanup()
})
