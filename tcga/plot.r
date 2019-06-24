library(dplyr)
library(patchwork)
sys = import('sys')
plt = import('plot')

do_plot = function(data, title) {
    df = data %>%
        mutate(label=gene, size=n_aneup) %>%
        plt$color$p_effect(pvalue="adj.p", effect="estimate")
    p1 = df %>%
        filter(estimate < 0) %>%
        plt$volcano(base.size=0.1, label_top=30, repel=TRUE) +
        labs(title = sprintf("%s (%i aneup genes)", title, sum(!is.na(df$estimate))),
             subtitle = "compensated")
    p2 = df %>%
        filter(estimate > 0) %>%
        plt$volcano(base.size=0.1, label_top=30, repel=TRUE) +
        labs(subtitle="hyper-deregulated") +
        theme(axis.title.y = element_blank())
    p1 + p2 + plot_layout(nrow=1)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'xlsx', 'naive/LUAD.xlsx'),
        opt('p', 'plotfile', 'pdf', 'naive/LUAD.pdf'))

    fits = readxl::excel_sheets(args$infile) %>%
        sapply(function(s) readxl::read_xlsx(args$infile, sheet=s), simplify=FALSE)
#    fits2 = unlist(lapply(fits, function(x) split(x, x$term)), recursive=FALSE)

    pdf(args$plotfile, 10, 8)
    for (i in seq_along(fits))
        print(do_plot(fits[[i]], names(fits)[i]))
    dev.off()
})
