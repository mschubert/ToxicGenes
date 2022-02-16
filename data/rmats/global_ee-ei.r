library(dplyr)
library(GenomicFeatures)
library(GenomicAlignments)
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'global_ee-ei.rds'),
    opt('p', 'plotfile', 'pdf', 'global_ee-ei.pdf')
)

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
meta$reads = lapply(reads, . %>% as.data.frame() %>% tibble::rownames_to_column("gene_name") %>% as_tibble())
saveRDS(meta, file=args$outfile)

mu = tidyr::unnest(meta, reads) #%>% filter(cline == "HCC70", time == "24h")

do_glm = function(df) {
#    glm(cbind(ee, ei) ~ cline + cond, data=df, family="binomial") %>%
    MASS::glm.nb(ei ~ cline:ee + cond:ee, df) %>%
        broom::tidy() %>%
#        filter(term == "condRBM14") %>%
        filter(term == "ee:condRBM14") %>%
        select(-term)
}
xx = mu %>%
    mutate(ee_ei = ee / ei) %>%
    group_by(time, gene_name) %>%
        tidyr::nest() %>%
    ungroup() %>%
    mutate(res = clustermq::Q(do_glm, df=data, n_jobs=10, pkgs="dplyr", fail_on_error=FALSE)) %>%
    filter(!sapply(res, function(r) identical(class(r), "error"))) %>%
    tidyr::unnest(res) %>%
#    filter(std.error < 5) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value, -abs(statistic))

gset = import('genesets')
plt = import('plot')
gset$get_human("GO_Biological_Process_2021") %>% gset$test_lm(xx %>% filter(time == "8h"), .) %>% plt$volcano()


mu2 = mu %>% filter(ee > 50, ee < 1000) %>% mutate(ee_ei = ee/ei)
mu2.8 = mu2 %>% filter(time == "8h", gene_name %in% sample(unique(gene_name), 1000))
mu2.24 =mu2 %>% filter(time == "24h", gene_name %in% sample(unique(gene_name), 1000))

mod8 = lm(ee_ei ~ gene_name + cond, data=mu2.8)
mod24 = lm(ee_ei ~ gene_name + cond, data=mu2.24)

broom::tidy(mod8) %>% filter(!grepl("gene_name", term))
broom::tidy(mod24) %>% filter(!grepl("gene_name", term))

lm(ee_ei ~ gene_name + cond, data=mu2.8 %>% filter(cline == "H838")) %>% broom::tidy() %>% filter(!grepl("gene_name", term))
lm(ee_ei ~ gene_name + cond, data=mu2.8 %>% filter(cline == "H1650")) %>% broom::tidy() %>% filter(!grepl("gene_name", term))
lm(ee_ei ~ gene_name + cond, data=mu2.8 %>% filter(cline == "HCC70")) %>% broom::tidy() %>% filter(!grepl("gene_name", term))

lm(ee_ei ~ gene_name + cond, data=mu2.24 %>% filter(cline == "H838")) %>% broom::tidy() %>% filter(!grepl("gene_name", term))
lm(ee_ei ~ gene_name + cond, data=mu2.24 %>% filter(cline == "H1650")) %>% broom::tidy() %>% filter(!grepl("gene_name", term))
lm(ee_ei ~ gene_name + cond, data=mu2.24 %>% filter(cline == "HCC70")) %>% broom::tidy() %>% filter(!grepl("gene_name", term))

pdf(args$plotfile, 10, 8)
ggplot(mu %>% filter(ee > 50), aes(x=rep, y=log2(ei/ee))) +
    geom_hline(yintercept=0, linetype="dashed", color="grey", size=1.5) +
    geom_boxplot(aes(fill=cond), outlier.shape=NA, alpha=0.6) +
    scale_fill_manual(values=c(Luc="black", RBM14="red")) +
    facet_grid(cline ~ time) +
    coord_cartesian(ylim=c(-3,3)) +
    theme_minimal()

ggplot(mu %>% filter(ee > 5, ei > 5), aes(x=ee, y=ei)) +
    geom_density2d(aes(color=cond, linetype=rep), bins=5) +
    scale_color_manual(values=c(Luc="black", RBM14="red")) +
    scale_linetype_manual(values=c("1"="dashed", "2"="dotted")) +
    facet_grid(cline ~ time) +
    scale_x_continuous(trans="log1p", breaks=10^(1:4), limits=c(5,1e4), expand=c(0,0)) +
    scale_y_continuous(trans="log1p", breaks=10^(1:4), limits=c(5,1e4), expand=c(0,0)) +
    theme_minimal()
dev.off()
