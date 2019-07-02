library(dplyr)
library(patchwork)
plt = import('plot')
sys = import('sys')

do_plot = function(data, title) {
    ylab = "Adjusted p-value (FDR)"
    if (sum(data$adj.p < 1e-300) > 5) {
        data$adj.p = 10^-abs(pmin(data$statistic, 300))
        ylab = "Pseudo p-value (values too close to zero)"
    }

    p1 = tryCatch({ 
        p = data %>%
            filter(estimate < 0) %>%
            plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
            plt$volcano(base.size=0.2, text.size=2.5, label_top=30, repel=TRUE) +
            labs(title = title,
                 subtitle = "compensated",
                 y = ylab)
        ggplot_build(p)
        p
    }, error = function(e) plt$text(conditionMessage(e)))
        
    p2 = tryCatch({
        p = data %>%
            filter(estimate > 0) %>%
            plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
            plt$volcano(base.size=0.2, text.size=2.5, label_top=30, repel=TRUE) +
            labs(subtitle="hyper-deregulated") +
            theme(axis.title.y = element_blank())
        ggplot_build(p)
        p
    }, error = function(e) plt$text(conditionMessage(e)))

    p1 + p2 + plot_layout(nrow=1)
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

    pdf(args$plotfile, 12, 8)
    for (i in seq_along(fits)) {
        message(names(fits)[i])
        print(do_plot(fits[[i]], names(fits)[i]))
    }
    dev.off()
})
