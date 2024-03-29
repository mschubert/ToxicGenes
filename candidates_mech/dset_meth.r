library(dplyr)
library(GenomicRanges)
sys = import('sys')
seq = import('seq')
plt = import('plot')
tcga = import('data/tcga')
util = import('../candidates/util')
util2 = import('./util')

get_ensembl_rest_info = function(gene_id) {
    query = sprintf("https://rest.ensembl.org/lookup/id/%s?expand=1", gene_id)
    r = httr::GET(query, httr::content_type("application/json"))
    httr::stop_for_status(r)
    res = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))

    res$Transcript %>% select(ensembl_transcript_id=id, Exon) %>% tidyr::unnest(Exon) %>%
        lapply(unlist) %>% do.call(data.frame, .) %>% as_tibble()
}

get_ensembl_regulatory = function(.chr, .start, .stop) {
    query = sprintf("https://rest.ensembl.org/overlap/region/human/%s:%i-%i?feature=regulatory", .chr, .start, .stop)
    r = httr::GET(query, httr::content_type("application/json"))
    httr::stop_for_status(r)
    res = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))

    res %>% select(label=id, start, end, text=description) %>%
        lapply(unlist) %>% do.call(data.frame, .) %>% as_tibble() %>%
        mutate(start = pmax(start, .start),
               end = pmin(end, .stop))
}

get_ensembl_motif = function(.chr, .start, .stop) {
    query = sprintf("https://rest.ensembl.org/overlap/region/human/%s:%i-%i?feature=motif", .chr, .start, .stop)
    r = httr::GET(query, httr::content_type("application/json"))
    httr::stop_for_status(r)
    res = jsonlite::fromJSON(jsonlite::toJSON(httr::content(r)))

    res %>% mutate(ev = ifelse(sapply(epigenomes_with_experimental_evidence, is.null),
                               NA, epigenomes_with_experimental_evidence)) %>%
        select(start, end, tf=transcription_factor_complex, score, ev) %>%
        lapply(unlist) %>% do.call(data.frame, .) %>% as_tibble() %>%
        mutate(tf = strsplit(tf, ",|::")) %>%
        tidyr::unnest("tf") %>%
        group_by(tf) %>% slice_max(score, n=3) %>% ungroup()
}

plot_gene_annot = function(gene_name, cpg) {
    trs = seq$coords$transcript() %>%
        filter(external_gene_name == gene_name)
    strand = setNames(c("-", "+"), c("-1","1"))[as.character(trs$strand[1])]
    adds = get_ensembl_rest_info(trs$ensembl_gene_id[1]) %>%
        select(label=ensembl_transcript_id, exon_id=id, start, end) %>%
        mutate(type="transcript")
    trs2 = trs %>%
        transmute(label=ensembl_transcript_id, start=transcript_start, end=transcript_end)

    depmap_tss = readr::read_tsv("../data/ccle/CCLE_RRBS_TSS_1kb_20180614.txt") %>%
        filter(gene == gene_name) %>%
        transmute(label=sprintf("DepMap_TSS_%i", seq_len(nrow(.))), chr=chr, start=fpos, end=tpos)
    depmap_cpg = readr::read_tsv("../data/ccle/CCLE_RRBS_TSS_CpG_clusters_20180614.txt") %>%
        dplyr::rename(gene = gene_name) %>%
        filter(gene == gene_name) %>%
        rowwise() %>%
        transmute(chr = paste0("chr", sub("([0-9]+):.*", "\\1", CpG_sites_hg19)),
                  start = strsplit(CpG_sites_hg19, ";")[[1]] %>% sub("[0-9]+:", "", .) %>% as.numeric() %>% min(),
                  end = strsplit(CpG_sites_hg19, ";")[[1]] %>% sub("[0-9]+:", "", .) %>% as.numeric() %>% max()) %>%
        ungroup() %>%
        mutate(label = sprintf("DepMap_CpG_%i", seq_len(nrow(.))))
    depmap = bind_rows(depmap_tss, depmap_cpg)
    gr = GenomicRanges::makeGRangesFromDataFrame(depmap, keep.extra.columns=TRUE)
    seqlevelsStyle(gr) = "UCSC"
    hg19_to_hg38 = rtracklayer::import.chain("../data/ccle/hg19ToHg38.over.chain")
    depmap = rtracklayer::liftOver(gr, hg19_to_hg38) %>% unlist() %>% as.data.frame() %>%
        mutate(cg_sd = 0.01)

    probes = grep("^cg", colnames(cpg), value=TRUE)
    cgs = tcga$meth_cg2gene() %>%
        filter(probeID %in% probes) %>%
        transmute(label=probeID, start=CpG_beg, end=CpG_end) %>%
        left_join(tibble(label=colnames(cpg), cg_sd=apply(cpg, 2, sd)))

    islands = tcga$meth_cg2island() %>% filter(probeID %in% probes) %>%
        mutate(changed = cumsum(CGIposition != c(CGIposition[1], CGIposition[-nrow(.)]))) %>%
        filter(CGIposition == "Island") %>%
        group_by(changed) %>%
        summarize(start=min(CpG_beg), end=max(CpG_end))

    both = bind_rows(list(transcript=trs2, cgs=cgs, DepMap=depmap), .id="type")
    .chr = unique(trs$chromosome_name)
    .start = min(both$start, na.rm=TRUE) - 1000
    .stop = max(both$end, na.rm=TRUE) + 1000
    reg = get_ensembl_regulatory(.chr, .start, .stop) %>% mutate(type="regulatory")
    tfs = get_ensembl_motif(.chr, .start, .stop)

    ccle_track = readr::read_tsv("../data/ccle/methylation_gene_body.tsv") %>%
        filter(chr == paste0("chr", .chr), position >= .start, position <= .stop) %>%
        group_by(position) %>%
        summarize(beta = mean(beta), n_clines=n_distinct(cell_line)) %>%
        transmute(label="CCLE David Wu", start=position, type="BS-seq", cg_sd=n_clines/(max(n_clines)*100))
    mdamb_track = rtracklayer::import(rtracklayer::BigWigFile("merged_file_sample.bw")) %>%
        as.data.frame() %>% as_tibble() %>%
        filter(seqnames %in% .chr, start >= .start, end <= .stop) %>%
        transmute(label="MDA-MB BS-seq", start=start, meth=score, type="BS-seq")
    both = bind_rows(list(both, reg, mdamb_track, ccle_track))

    ggplot(both, aes(x=start, y=label, color=type)) +
        geom_rect(data=islands, aes(xmin=start, xmax=end), ymin=-Inf, ymax=Inf,
                  inherit.aes=FALSE, fill="black", alpha=0.1) +
        ggridges::geom_ridgeline(aes(height=meth/50, fill=type)) +
        geom_segment(aes(xend=end, yend=label), size=1) +
        geom_segment(data=adds, aes(xend=end, yend=label), size=4, alpha=0.5) +
        geom_point(aes(size=cg_sd)) +
        geom_text(aes(label=text), size=3, color="black") +
        scale_size_area(max_size=8) +
        scale_x_continuous(expand=c(0,0)) +
        ggrepel::geom_text_repel(data=tfs, aes(y=rank(score)/length(score)*length(unique(both$label)),
                                               label=tf, fontface=ifelse(is.na(ev), "plain", "bold")),
            color="black", size=ifelse(is.na(tfs$ev), 1, 2.5), alpha=0.7, segment.alpha=0.1, max.overlaps=Inf) +
        ggtitle(sprintf("%s (%s)", gene_name, strand))
}

meth_for_gene = function(cohorts, gene_name, gene_id="external_gene_name") {
    cg1 = tcga$meth_summary(cohorts, gene_id)[gene_name,,]
    map = tcga$meth_mapping(gene_id)$pgene
    probe_ids = map$probe_id[map[[gene_id]] == gene_name]
    cg2 = lapply(cohorts, function(c) {
        m = tcga$meth_cpg(c)
        m[intersect(rownames(m), probe_ids),,drop=FALSE]
    }) %>% narray::stack(along=2) %>% t()
    cg2 = cg2[,names(sort(-apply(cg2, 2, sd)))]
    re = narray::stack(cg1, cg2, along=2)
    re[,colSums(is.na(re)) < nrow(re)/2]
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'gene', 'HGNC symbol', 'CDKN1A'),
        opt('i', 'infile', 'rds', 'by_gene/CDKN1A.rds'),
        opt('p', 'plotfile', 'pdf', 'CDKN1A/meth.pdf')
    )

    td = readRDS(args$infile)
    cohorts = unique(td$cohort) %>% setdiff(c("BRCA.LumAB", "BRCA.Basal"))
    cpg = tryCatch(meth_for_gene(cohorts, args$gene), error=util2$muffle)

    td = td %>% filter(sample %in% rownames(cpg)) %>%
        mutate(cohort = ifelse(cohort %in% c("COAD", "READ"), "COADREAD", cohort),
               cohort = ifelse(cohort %in% c("LUAD", "LUSC"), "NSCLC", cohort))
    cpg = cpg[td$sample,]
    dset = cbind(td, cpg, constant=1)

    pdf(args$plotfile, 28, 12)
    print(plt$text(sprintf("Methylation (%i CpG)", util2$nc(cpg)), size=20))
    print(plot_gene_annot(args$gene, cpg))
    if (sum(!is.na(dset$meth_eup_scaled)) > 0)
        print(util2$plot_l2d(dset, "meth_eup_scaled"))

    for (v in colnames(cpg))
        print(util2$plot_l2d(dset, v))
    dev.off()
})
