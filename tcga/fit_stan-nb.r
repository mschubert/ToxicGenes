library(dplyr)
sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

#' Fit a negative binomial model using a stan glm
#'
#' @param df    A data.frame with columns: purity, expr, copies, sf, covar, eup_equiv
#' @param cna   Character string of either: "amp", "del", or "all"
#' @param type  Character string of either: "naive", "pur"
#' @param et    Tolerance in ploidies to consider sample euploid
do_fit = function(df, cna, type="pur", et=0.15) {
    if (cna == "amp") {
        df = df %>% filter(copies > 2-et)
    } else if (cna == "del") {
        df = df %>% filter(copies < 2+et)
    }

    if (length(unique(df$covar)) == 1) {
        full = switch(type,
            naive = expr ~ eup_equiv,
            pur = expr ~ purity + eup_equiv
#            puradj = expr ~ purity + eup_equiv
        )
    } else {
        full = switch(type,
            naive = expr ~ covar + eup_equiv,
            pur = expr ~ covar + purity + eup_equiv
#            puradj = expr ~ covar + purity + eup_equiv
        )
    }

    tryCatch({
        res = rstanarm::stan_glm(full, data=df, offset=sizeFactor,
                                 family = neg_binomial_2(link="identity"),
                                 prior = normal(mean(df$expr), sd(df$expr)),
                                 prior_intercept = normal(mean(df$expr), sd(df$expr)),
                                 seed = 2380719)

        rmat = as.matrix(res)
        hdist = statip::hellinger(rmat[,"(Intercept)"], rmat[,"eup_equiv"])

        pseudo_p = pnorm(abs((mean(rmat[,"(Intercept)"]) - mean(rmat[,"eup_equiv"])) /
                             sd(rmat[,"(Intercept)"])), lower.tail=F)

        tibble(estimate = mean(rmat[,"eup_equiv"]) / mean(rmat[,"(Intercept)"]) - 1, # pct_comp
               n_aneup = sum(abs(df$copies-2) > 1-et),
               n_genes = 1,
               eup_reads = mean(rmat[,"(Intercept)"]),
               slope_diff = mean(rmat[,"(Intercept)"]) - mean(rmat[,"eup_equiv"]),
               rsq = hdist, # not rsq, but [0,1]
               p.value = pseudo_p)

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
        opt('j', 'cores', 'integer', '10'),
        opt('m', 'memory', 'integer', '4096'),
        opt('o', 'outfile', 'xlsx', 'LUAD/genes.xlsx')
    )

    et = yaml::read_yaml(args$config)$euploid_tol

    if (args$tissue == "pan")
        args$tissue = tcga$cohorts()

    # excl x,y chroms

    purity = tcga$purity() %>%
        filter(!is.na(estimate)) %>%
        transmute(sample = Sample,
                  purity = estimate)

    reads = lapply(args$tissue, tcga$rna_seq) %>%
        narray::stack(along=2) %>%
        tcga$filter(cancer=TRUE, primary=TRUE)
    rownames(reads) = idmap$gene(rownames(reads), to="hgnc_symbol")
    reads = reads[rowMeans(reads) >= 10 & !is.na(rownames(reads)),]
    eset = DESeq2::DESeqDataSetFromMatrix(reads, tibble(sample=colnames(reads)), ~1) %>%
        DESeq2::estimateSizeFactors() #todo: 16 duplicate hgnc_symbol

    copies = lapply(args$tissue, tcga$cna_genes) %>%
        narray::stack(along=2)
    rownames(copies) = idmap$gene(rownames(copies), to="hgnc_symbol")
    narray::intersect(reads, copies, along=1)
    narray::intersect(reads, copies, along=2)
    names(dimnames(reads)) = names(dimnames(copies)) = c("gene", "sample")

    df = inner_join(reshape2::melt(reads, value.name="expr"),
                    reshape2::melt(copies, value.name="copies")) %>%
        inner_join(purity) %>%
        inner_join(as.data.frame(SummarizedExperiment::colData(eset))) %>%
        na.omit() %>%
        mutate(covar = tcga$barcode2study(sample),
               cancer_copies = pmax(0, (copies-2) / purity + 2),
               eup_equiv = (cancer_copies - 2) / 2) %>%
        group_by(gene) %>%
        tidyr::nest() %>%
        tidyr::expand_grid(tibble(cna = c("amp", "del", "all")))

    rm(list=c("eset", "reads", "copies")); gc()

    df$res = clustermq::Q(do_fit, df=df$data, cna=df$cna,
                          const = list(type=args$type, et=et),
                          n_jobs = as.integer(args$cores), memory=args$memory,
                          pkgs = c("dplyr", "rstanarm"), chunk_size=1)

    res = df %>%
        select(-data) %>%
        tidyr::unnest("res") %>%
        split(.$cna) %>%
        lapply(. %>% mutate(adj.p = p.adjust(p.value, method="fdr")) %>% arrange(adj.p, p.value))

    writexl::write_xlsx(res, args$outfile)
})
