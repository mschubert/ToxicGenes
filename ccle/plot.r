library(dplyr)
plt = import('plot')
sys = import('sys')

do_plot = function(data) {
    data %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
        plt$volcano(base.size=0.2, label_top=50, repel=TRUE,
                    text.size=2.5, x_label_bias=5, pos_label_bias=0.15)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'xlsx', 'pan/genes.xlsx'),
        opt('p', 'plotfile', 'pdf', 'pan.pdf'))

    fits = readxl::excel_sheets(args$infile) %>%
        sapply(function(s) readxl::read_xlsx(args$infile, sheet=s), simplify=FALSE)

    if (grepl("genes\\.xlsx", args$infile))
        fits = lapply(fits, function(f) mutate(f, label=gene, size=n_aneup))
    else
        fits = lapply(fits, function(f) mutate(f, label=set))

    pdf(args$plotfile)
    for (i in seq_along(fits))
        print(do_plot(fits[[i]]) + ggtitle(names(fits)[i]))
    dev.off()
})
