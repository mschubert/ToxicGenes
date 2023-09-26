library(dplyr)
library(patchwork)
sys = import('sys')
gset = import('data/genesets')
plt = import('plot')

do_plot = function(res, sets, title="") {
    adjlab = function(lab, sz) {
        name_has_genes = all(sapply(sets[[lab]], function(s) grepl(s, lab, fixed=TRUE)))
        if (sz < 10 & ! name_has_genes)
            sprintf("%s\n%s", lab, paste(sets[[lab]], collapse=","))
        else
            sprintf("%s (%i)", lab, sz)
    }
    res = rowwise(res) %>%
        mutate(label = adjlab(label, size_used)) %>%
        ungroup()

    left = res %>%
        filter(estimate < 0) %>%
        plt$volcano(base.size=0.2, text.size=2.5, p=0.15, label_top=15, repel=TRUE)
    right = res %>%
        filter(estimate > 0) %>%
        plt$volcano(base.size=0.2, text.size=2.5, p=0.15, label_top=15, repel=TRUE)
    plt$text(title, size=6) / (left | right) + plot_layout(heights=c(1,20))
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('i', 'infile', 'xlsx', 'pan/stan-nb.xlsx'),
        opt('s', 'setfile', 'rds', '../data/genesets/CORUM_core.rds'),
        opt('e', 'select', 'all|comp|comp+orf', 'all'),
        opt('o', 'outfile', 'rds', 'pan/stan-nb/MSigDB_Hallmark_2020.rds'),
        opt('p', 'plotfile', 'pdf', 'pan/stan-nb/MSigDB_Hallmark_2020.pdf')
    )

    cnas = yaml::read_yaml(args$config)$cna

    dset = readxl::excel_sheets(args$infile) %>%
        sapply(function(s) readxl::read_xlsx(args$infile, sheet=s), simplify=FALSE)

    sets = readRDS(args$setfile)

    sel = readr::read_tsv("../cor_tcga_ccle/positive_comp_set.tsv") %>% filter(hit)
    if (args$select == "comp+orf")
        sel = sel %>% filter(est_orf < -1)
    if (grepl("comp", args$select))
        sets = sets[sapply(sets, function(s) any(sel$gene %in% s))]

    result = lapply(dset, function(d) gset$test_lm(d, sets, stat="estimate"))
    plots = mapply(do_plot, res=result, title=names(result),
                   MoreArgs=list(sets=sets), SIMPLIFY=FALSE)

    pdf(args$plotfile, 12, 8)
    for (p in plots)
        print(p)
    dev.off()

    saveRDS(result, file=args$outfile)
})
