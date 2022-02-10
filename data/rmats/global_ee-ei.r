library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)

genomeTxDb = makeTxDbFromGFF("seqdata/hg38.refseq.gtf")
bams = list.files("seqdata", "\\.bam$", recursive=TRUE, full.names=TRUE)

quant_one = function(fname) {
    ga = GenomicAlignments::readGAlignments(Rsamtools::BamFile(fname), use.names=TRUE,
                index=paste0(fname, ".bai"),
                param=Rsamtools::ScanBamParam(flag=Rsamtools::scanBamFlag(isDuplicate=FALSE)))

    juncs = as_tibble(cbind(#n = njunc(ga), cigarOpTable(cigar(ga)),
                            ee_junc = grepl("[0-9]+M[0-9]+N[0-9]+M", cigar(ga)),
                            ei_junc = grepl("([0-9]+M[0-9]+S)|([0-9]+S[0-9]+M)", cigar(ga))))

    cbind(ee = countOverlaps(genes(genomeTxDb), ga[juncs$ee_junc]),
          ei = countOverlaps(genes(genomeTxDb), ga[juncs$ei_junc]))
}

reads = parallel::mclapply(bams, quant_one, mc.cores=min(5, length(bams)))
#names(reads) = sub("\\.sorted\\.bam", "", basename(h838))

meta = sub(".sorted.bam", "", basename(bams), fixed=TRUE) %>%
    strsplit("_") %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    transmute(cline = factor(V1),
              cond = sub("luc", "Luc", V2) %>% factor() %>% relevel("Luc"),
              time = factor(V3) %>% relevel("8h"),
              rep = factor(V4)) %>%
    as_tibble()
meta$reads = lapply(reads, as_tibble)

mu = tidyr::unnest(meta) #%>% filter(cline == "HCC70", time == "24h")

ggplot(mu %>% filter(ee > 1000), aes(x=rep, y=ei/ee)) +
    geom_boxplot(aes(fill=cond), outlier.shape=NA) +
    facet_grid(cline ~ time) +
    coord_cartesian(ylim=c(0,2))

ggplot(mu %>% filter(ee > 5, ei > 5), aes(x=ee, y=ei)) +
#    geom_point(aes(color=cond, shape=rep), alpha=0.02) +
    geom_density2d(aes(color=cond, linetype=rep), bins=5) +
    scale_color_manual(values=c(Luc="grey", RBM14="red")) +
    facet_grid(cline ~ time) +
    scale_x_continuous(trans="log1p", breaks=10^(1:4), limits=c(5,1e4), expand=c(0,0)) +
    scale_y_continuous(trans="log1p", breaks=10^(1:4), limits=c(5,1e4), expand=c(0,0))
