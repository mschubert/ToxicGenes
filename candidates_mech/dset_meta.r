library(dplyr)
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')
util2 = import('./util')

surv_knn = function(expr, cancer_copies, os_status, os_days, k=20) {
    death50 = function(ii) {
        cur = na.omit(os[ii,]) %>% arrange(days)
        n_deaths = cumsum(cur$status == "dead")
        cur$days[which(n_deaths >= rev(n_deaths)[1]/2)[1]]
    }
    mat = cbind(expr, cancer_copies)
    os = data.frame(status=os_status, days=os_days)
    noNA = !is.na(mat[,1]) & !is.na(mat[,2])
    kmat = FNN::get.knn(mat[noNA,], k=k)$nn.index

    re = rep(NA, nrow(mat))
    re[noNA] = apply(cbind(seq_len(nrow(kmat)), kmat), 1, death50)
    pmin(re, 365*5)
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/meta.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort) %>% setdiff(c("BRCA.LumAB", "BRCA.Basal"))

    aneup1 = lapply(cohorts, tcga$aneuploidy) %>%
        bind_rows() %>%
        select(Sample, aneuploidy)
    aneup2 = tcga$cna_absolute(cohorts) %>%
        mutate(cohort = tcga$barcode2study(Sample)) %>%
        group_by(cohort, Sample) %>%
            summarize(abs_ploidy = weighted.mean(Modal_HSCN_1 + Modal_HSCN_2, Length),
                      abs_aneup = weighted.mean(abs(Modal_HSCN_1 + Modal_HSCN_2 - 2), Length)) %>%
        ungroup()
    aneup = full_join(aneup1, aneup2) %>% dplyr::rename(sample=Sample)

    dset = td %>%
        mutate(constant = 1) %>%
        group_by(cohort, p53_mut) %>%
            mutate(death50_k5 = surv_knn(expr, cancer_copies, os_status, os_days, k=5),
                   death50_k20 = surv_knn(expr, cancer_copies, os_status, os_days, k=20)) %>%
        ungroup() %>%
        left_join(aneup)

    pdf(args$plotfile, 28, 12)
    print(util2$plot_l2d(dset, "purity", from=0.5, to=1, by="constant"))
    print(util2$plot_l2d(dset, "expr", from=0, by="constant"))
    print(util2$plot_l2d(dset, "death50_k5", by="constant"))
    print(util2$plot_l2d(dset, "death50_k20", by="constant"))
    print(util2$plot_l2d(dset, "aneuploidy", by="purity", RdBu=TRUE))
    print(util2$plot_l2d(dset, "abs_ploidy", to=4, by="purity", RdBu=TRUE))
    print(util2$plot_l2d(dset, "abs_aneup", to=2, by="purity", RdBu=TRUE))
    dev.off()
})
