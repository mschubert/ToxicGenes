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
            plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
            plt$volcano(base.size=0.1, text.size=2.5, label_top=30, repel=TRUE) +
            labs(title = title,
                 subtitle = "compensated")
        plt$try(p)
    }, error = function(e) plt$text(conditionMessage(e)))
        
    p2 = tryCatch({
        p = data %>%
            filter(estimate > 0) %>%
            plt$color$p_effect(pvalue="adj.p", effect="estimate") %>%
            plt$volcano(base.size=0.1, text.size=2.5, label_top=30, repel=TRUE) +
            labs(subtitle="hyper-deregulated") +
            theme(axis.title.y = element_blank())
        plt$try(p)
    }, error = function(e) plt$text(conditionMessage(e)))

    p1 + p2 + plot_layout(nrow=1)
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'xlsx', 'pan/rlm3.xlsx'),
        opt('p', 'plotfile', 'pdf', 'pan.pdf'))

    fits = readxl::excel_sheets(args$infile) %>%
        sapply(function(s) readxl::read_xlsx(args$infile, sheet=s), simplify=FALSE)
    fits = lapply(fits, function(f) mutate(f, label=name, size=n_aneup))

    pdf(args$plotfile, 12, 8)
    for (i in seq_along(fits)) {
        message(names(fits)[i])
        print(do_plot(fits[[i]], names(fits)[i]))
    }
    dev.off()
})
