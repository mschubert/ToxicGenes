library(ggplot2)
library(dplyr)
library(plyranges)
library(Gviz)
options(ucscChromosomeNames=FALSE)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('g', 'gtf', 'assembly GTF', 'seqdata/hg38.refseq.gtf'),
    opt('f', 'flank', 'integer', '1000'),
    opt('p', 'plotfile', 'PDF to save to', 'sashimi.pdf')
)

bams = list.files("seqdata", "\\.bam$", recursive=TRUE, full.names=TRUE)
meta = sub(".sorted.bam", "", basename(bams), fixed=TRUE) %>%
    strsplit("_") %>%
    do.call(rbind, .) %>%
    as.data.frame() %>%
    transmute(cline = factor(V1),
              cond = sub("luc", "Luc", V2) %>% factor() %>% relevel("Luc"),
              time = factor(V3) %>% relevel("8h"),
              rep = factor(V4),
              bam = bams) %>%
    as_tibble() %>%
    filter(time != "72h")

txdb = GenomicFeatures::makeTxDbFromGFF(args$gtf, format="gtf")
genes = seq$coords$gene(assembly="GRCh38", granges=TRUE)
exons = GenomicFeatures::exonsBy(txdb, by="tx", use.names=TRUE)

.name = "NEK7"
.time = "24h"

region = genes[genes$external_gene_name %in% .name] %>%
    anchor_5p() %>% stretch(as.integer(args$flank)) %>%
    anchor_3p() %>% stretch(as.integer(args$flank))
seqlevelsStyle(region) = "UCSC"

grtrack = GeneRegionTrack(txdb, name=as.character(seqnames(region)))
gtrack = GenomeAxisTrack(name="GRCh38_RefSeq")
tracks = list(gtrack, grtrack)

rbm = tracks
cur = meta %>% filter(cond == "RBM14", time == .time) %>%
    mutate(name=paste(cline, cond, time, rep, sep="_")) %>%
    rowwise() %>%
        mutate(rna = list(AlignmentsTrack(bam, type=c("coverage","sashimi"), name=name))) %>%
    ungroup()
for (i in seq_len(nrow(cur)))
    rbm[[cur$name[i]]] = cur$rna[[i]]

luc = tracks
cur = meta %>% filter(cond == "Luc", time == .time) %>%
    mutate(name=paste(cline, cond, time, rep, sep="_")) %>%
    rowwise() %>%
        mutate(rna = list(AlignmentsTrack(bam, type=c("coverage","sashimi"), name=name))) %>%
    ungroup()
for (i in seq_len(nrow(cur)))
    luc[[cur$name[i]]] = cur$rna[[i]]
sizes = c(1,1,rep(2,nrow(cur)))

pdf("test.pdf", 6, 8)
plotTracks(rbm, from=start(region), to=end(region), chromosome=seqnames(region),
           sizes=sizes, cex.main=1, main=paste(.name, "RBM14", .time))
plotTracks(luc, from=start(region), to=end(region), chromosome=seqnames(region),
           sizes=sizes, cex.main=1, main=paste(.name, "Luc", .time))
dev.off()



eset = readRDS("eset_exon.rds") %>% DESeq2::estimateSizeFactors() %>% DESeq2::counts(normalized=TRUE)
expr = eset[grepl("^NEK7", rownames(eset)),] %>%
    reshape2::melt() %>%
    transmute(exon=Var1, bam=Var2, counts=value) %>%
    inner_join(meta %>% mutate(bam = basename(bam))) %>%
    as_tibble() %>%
    mutate(bam = sub("\\.sorted\\.bam", "", bam))

ggplot(expr, aes(x=exon, y=bam, fill=log2(counts+1))) +
    geom_raster() +
    facet_grid(cline ~ ., scales="free_y", space="free_y") +
    scale_fill_brewer(pal="")
