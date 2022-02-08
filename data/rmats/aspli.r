library(dplyr)
library(ASpli)
library(GenomicFeatures)
library(plyranges)

genomeTxDb = makeTxDbFromGFF("seqdata/hg38.refseq.gtf")
features = binGenome(genomeTxDb)

h838 = list.files("seqdata", "H838.*8h.*\\.bam$", recursive=TRUE, full.names=TRUE)
targets = data.frame(bam=h838, f1=ifelse(grepl("RBM14", h838), "treatment", "control"))

gbcounts = gbCounts(features=features, targets=targets, minReadLength=100, maxISize=50000)
asd = jCounts(counts=gbcounts, features=features, minReadLength=100)




# todo: read bam file, count e-e, e-i junctions only per gene
seq = import('seq')
ex = exonsBy(genomeTxDb, by="gene") %>% unlist() %>% mutate(gene_name = names(.)) %>% unname()
fname = h838[1]

process_one = function(fname) {
    gr = seq$read_granges(fname) %>%
        mutate(read_name = names(.)) %>% unname()

    within_ex = join_overlap_inner_within(gr, ex)
    spanning_ex = join_overlap_inner(gr, ex) %>%
        group_by(read_name) %>%
            filter(n_distinct(gene_name) == 1) %>%
        ungroup()



    juncs = join_overlap_left(gr, ex) %>%
        filter(!is.na(exon_id)) %>%
        filter(! read_name %in% within_ex$read_name)

    ee = juncs %>%
        group_by(gene_name, read_name) %>%
            filter(n_distinct(exon_id) > 1) %>%
        ungroup()

    ei = juncs %>%
        filter(! read_name %in% ee$read_name)

    list(ee=ee, ei=ei)
}
