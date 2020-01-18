library(ggplot2)
library(patchwork)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('t', 'tissue', 'pan|TCGA identifier', 'pan'),
    opt('f', 'fit', 'fit type id', 'rlm3'),
    opt('b', 'briscut', 'txt file', '../data/all_BrISCUT_results.txt'),
    opt('r', 'ranks', 'xlsx', 'rank_top/pan_rlm3.xlsx'),
    opt('p', 'plotfile', 'pdf', 'compare_genomic/pan_rlm3.pdf'))

top = sapply(readxl::excel_sheets(args$ranks), simplify=FALSE,
             function(s) readxl::read_xlsx(args$ranks, sheet=s)) %>%
    bind_rows(.id="cnv")

cnv = setNames(c(1, -1), c("amp", "del"))
genomic = readr::read_tsv(args$briscut) %>%
    transmute(cohort=type, name=Gene, log10_ks_p=pmin(log10_ks_p, 10) * cnv[direction])
if (args$tissue != "pan")
    genomic = genomic %>% filter(cohort == args$tissue)
genomic = genomic %>% group_by(name) %>% summarize(log10_ks_p = mean(log10_ks_p))

dset = inner_join(top, genomic, by="name") %>%
    group_by(cnv) %>%
    mutate(label = ifelse(rank(-score) < 30, name, NA),
           alpha = ifelse(is.na(label), 0.1, 0.8)) %>%
    ungroup()

pdf(args$plotfile, 8, 12)

ggplot(dset, aes(x=log10_ks_p, y=score)) +
    geom_point(aes(alpha=alpha)) +
    geom_vline(xintercept=0, linetype="dashed", size=1, color="darkgreen") +
    ggrepel::geom_text_repel(aes(label=label)) +
    facet_wrap(~ cnv, ncol=1) +
    guides(alpha=FALSE)

dev.off()
