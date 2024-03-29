library(dplyr)
library(patchwork)
sys = import('sys')
plt = import('plot')
deseq = import('process/deseq')
gset = import('genesets')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'depmap.rds'),
    opt('p', 'plotfile', 'pdf', 'depmap_sets.pdf')
)

sets = c("MSigDB_Hallmark_2020", "DoRothEA", "GO_Biological_Process_Tree", "CORUM_all")

idx = readRDS(args$infile) %>%
    filter(dset %in% c("rnai", "crispr_ko")) %>%
    dplyr::rename(genes = res) %>%
    deseq$sets(sets, cl=5)

hl = gset$get_human("MSigDB_Hallmark_2020")[["Oxidative Phosphorylation"]]
idx$genes = lapply(idx$genes, . %>% mutate(circle = label %in% hl))

plots = idx %>%
    tidyr::pivot_longer(all_of(c("genes", sets))) %>%
    rowwise() %>%
    mutate(plot = list(plt$volcano(value) + ggtitle(sprintf("%s (%s)", name, cond)))) %>%
    group_by(dset, field) %>%
    summarize(asm = list((plt$text(sprintf("%s :: %s", dset[1], field[1])) /
                          wrap_plots(plot, nrow=2)) + plot_layout(heights=c(1,20))))

pdf(args$plotfile, 22, 10)
for (p in plots$asm)
    print(p)
dev.off()
