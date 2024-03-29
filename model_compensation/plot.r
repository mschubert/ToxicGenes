library(dplyr)
library(patchwork)
plt = import('plot')
sys = import('sys')

do_plot = function(data, title) {
    fps = !is.na(data$estimate) & data$estimate < -5
    data$name[fps] = sprintf("%s:%i", data$name, as.integer(data$estimate))[fps]
    data$estimate[fps] = -5

    fps = !is.na(data$estimate) & data$estimate > 10
    data$name[fps] = sprintf("%s:%i", data$name, as.integer(data$estimate))[fps]
    data$estimate[fps] = 10

    p1 = tryCatch({
        p = data %>%
            filter(estimate < 0) %>%
            plt$volcano(base.size=0.1, text.size=2.5, label_top=30, repel=TRUE) +
                labs(title = title,
                     subtitle = "compensated")
        plt$try(p)
    }, error = function(e) plt$text(conditionMessage(e)))

    p2 = tryCatch({
        p = data %>%
            filter(estimate > 0) %>%
            plt$volcano(base.size=0.1, text.size=2.5, label_top=30, repel=TRUE) +
                labs(subtitle="hyper-deregulated") +
                theme(axis.title.y = element_blank())
        plt$try(p)
    }, error = function(e) plt$text(conditionMessage(e)))

    p1 + p2 + plot_layout(nrow=1)
}

rsq_vs_comp = function(data, title) {
    d2 = data %>%
        filter(slope_diff < 0, rsq > 0) %>%
        mutate(rnk2 = rank(slope_diff * rsq),
               rnk3 = rank(slope_diff),
               rnk4 = rank(-rsq),
               label = ifelse(rnk2 < 30 | rnk3 < 10 | rnk4 < 10,
                              sprintf("%s (%i)", name, n_aneup), NA))
    ggplot(d2, aes(x=-slope_diff, y=rsq, size=n_aneup)) +
        geom_point(alpha=0.1) +
        geom_hline(yintercept=0.05, linetype="dashed") +
        geom_vline(xintercept=0.5, linetype="dashed") +
        ggrepel::geom_label_repel(aes(label=label), size=3, label.size=NA,
            segment.alpha=0.3, min.segment.length=0, fill="#ffffffc0", label.padding=0.1) +
        theme_classic() +
        labs(title = title, subtitle="compensated")
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'xlsx', 'fit_brms/pan.rds'),
        opt('p', 'plotfile', 'pdf', 'plot_brms/pan.pdf'))

    fits = readRDS(args$infile)
    fits = lapply(fits, function(f) {
        if (! "name" %in% colnames(f)) # 2 lines: tcga meth 'gene' as colname
            f$name = f$gene
        mutate(f, label=name, size=n_aneup)
    })

    pdf(args$plotfile, 12, 8)
    for (i in seq_along(fits)) {
        message(names(fits)[i])
        print(do_plot(fits[[i]], names(fits)[i]))
        if ("rsq" %in% colnames(fits[[i]]))
            print(rsq_vs_comp(fits[[i]], names(fits)[i]))
    }
    dev.off()
})
