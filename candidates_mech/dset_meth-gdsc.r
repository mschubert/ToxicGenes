library(dplyr)
library(GenomicRanges)
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
gdsc = import('data/gdsc')
util = import('../candidates/util')
ccle_meth = import('./dset_meth-ccle')

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/meth-gdsc.pdf')
    )

    stopifnot(args$gene == "CDKN1A") # p21 hardcoded for now

    cohorts = c("BRCA", "COADREAD", "LUAD", "LUSC", "SKCM")
    cd = util$load_ccle(c(args$gene,"TP53")) %>%
        filter(gene == args$gene,
               cohort %in% cohorts) %>%
        select(-meth) %>%
        mutate(cosmic_id = gdsc$cosmic$name2id(Name)) %>%
        filter(!is.na(cosmic_id)) # join w/ cg explodes otherwise

    cgs = c("cg13662121", "cg11920449", "cg24425727", "cg15474579", "cg03032677", "cg01955533")
    cg = gdsc$meth()[cgs,] %>% #TODO:
        reshape2::melt() %>%
        as_tibble() %>%
        transmute(meth_id=Var1, meth_value=value, cosmic_id=as.character(Var2))

    comb = inner_join(cd, cg)

    pdf(args$plotfile, 16, 15)
    print(ccle_meth$plot_ccle_meth(comb))
    dev.off()
})
