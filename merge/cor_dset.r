library(dplyr)
library(ggplot2)
b = import('base')
sys = import('sys')

#' Plot compensation/hyperactivation between two data sets
#'
#' @param a1     Character string of data set on x axis
#' @param a2     Character string of data set on y axis
#' @param smat   Wald statistic matrix [genes x data sets]
#' @param cap    If values fall outside cap, plot them at the border
#' @param wald   Wald statistic to consider minimum for change (otherwise grey)
#' @param label  Number of points to label (consistent, inconsistent, orf, 1 ds)
#' @param cap_n  How many outliers must be outside cap in order to extend it
do_plot = function(a1, a2, smat, cap=20, wald=1.5, label=c(20, 10, 20, 3), cap_n=5) {
    message(a1, " & ", a2)
    x1 = rlang::sym(a1)
    x2 = rlang::sym(a2)
    colors = setNames(c("#a50f15", "#006d2c", "#045a8d", "#cccccc", "#8a8a8a", "#ffa81e"),
        c("compensated", "hyper-dereg", "inconsistent", "no change", "only 1 dset", "immune"))
    cur = smat[,c(a1, a2)]
    cap_at = function(x, m) sign(x) * pmin(abs(x), max(cap, -x[m][rank(x[m])==cap_n]))
    cur[,1] = cap_at(cur[,1], !is.na(cur[,2]))
    cur[,2] = cap_at(cur[,2], !is.na(cur[,1]))

    data = as.data.frame(cur) %>%
        tibble::rownames_to_column("name") %>%
        filter(! (is.na(!! x1) & is.na(!! x2)))

    plot_data = data %>%
        mutate(type = case_when(
                    grepl("TR[ABG]V|IG.V", name) ~ "immune",
                    is.na(!! x1) | is.na(!! x2) ~ "only 1 dset",
                    abs(!! x1) < wald | abs(!! x2) < wald ~ "no change",
                    !! x1 < 0 & !! x2 < 0 ~ "compensated",
                    !! x1 > 0 & !! x2 > 0 ~ "hyper-dereg",
                    TRUE ~ "inconsistent"
               ),
               type = factor(type, levels=names(colors)),
               topORF = name %in% rownames(smat)[rank(smat[,"orf"]) <= label[3]],
               fface = ifelse(topORF, "bold", "plain")) %>%
        group_by(sign(!! x1), sign(!! x2)) %>%
        mutate(annot_score = abs(!! x1 / quantile(!! x1, 0.9, na.rm=TRUE))^0.5 +
                             abs(!! x2 / quantile(!! x2, 0.9, na.rm=TRUE))^0.5,
               label = case_when(
                    type == "immune" ~ name,
                    topORF ~ name,
                    is.na(!! x1) & rank(-abs(!! x2)) <= label[4] ~ name,
                    is.na(!! x2) & rank(-abs(!! x1)) <= label[4] ~ name,
                    sign(!! x1) * sign(!! x2) > 0 & rank(-annot_score) <= label[1] ~ name,
                    sign(!! x1) * sign(!! x2) < 0 & rank(-annot_score) <= label[2] ~ name,
                    TRUE ~ as.character(NA)
               )) %>%
        ungroup()
    pcor = broom::tidy(lm(plot_data[[a1]] ~ plot_data[[a2]]))$p.value[2]

    nums = data %>%
        group_by(a1=sign(!! x1), a2=sign(!! x2)) %>%
        summarize(n = n()) %>%
        na.omit() %>%
        mutate(hjust=-a1/1.8+0.5, vjust=-a2+0.5)
    pfet = broom::tidy(fisher.test(matrix(pull(nums, n), ncol=2))) %catch%
        list(estimate=NA, p.value=NA)
    nums2 = data %>%
        group_by(a1=sign(!! x1), a2=sign(!! x2)) %>%
        filter(abs(!! x1) >= wald & abs(!! x2) >= wald) %>%
        summarize(n = n()) %>%
        na.omit() %>%
        mutate(x=sign(a1)*2, y=sign(a2)*2, hjust=-a1/2+0.5, vjust=-a2/2+0.5)
    pfet2 = broom::tidy(fisher.test(matrix(pull(nums2, n), ncol=2))) %catch%
        list(estimate=NA, p.value=NA)

    subt = sprintf(paste(sep = " | ",
                         "lm p=%.2g",
                         "comp. FET OR %.2f* (%.2f) p=%.2g* (%.2g)",
                         "stat cap %g",
                         "wald %g",
                         "top %i orf bold"),
                   pcor, pfet2$estimate, pfet$estimate,
                   pfet2$p.value, pfet$p.value, cap, wald, label[3])

    oneDS.x = min(plot_data[[a1]], na.rm=TRUE) - 1
    oneDS.y = min(plot_data[[a2]], na.rm=TRUE) - 1
    missing = filter(plot_data, is.na(!! x1) | is.na(!! x2))
    missing[[a1]][is.na(missing[[a1]])] = oneDS.x
    missing[[a2]][is.na(missing[[a2]])] = oneDS.y

    ggplot(plot_data, aes_string(x=a1, y=a2)) +
        geom_point(data=missing, aes(color=type)) +
        geom_point(aes(color=type), alpha=0.6) +
        geom_hline(yintercept=0, size=1, linetype="dashed", color="#dedede") +
        geom_vline(xintercept=0, size=1, linetype="dashed", color="#dedede") +
        scale_color_manual(values=colors) +
        geom_smooth(method="lm", color="black", linetype="dotted", se=FALSE, na.rm=TRUE) +
        geom_text(data=nums, aes(hjust=hjust, vjust=vjust, label=n), x=0, y=0, size=3) +
        geom_text(data=nums2, aes(x=x, y=y, hjust=hjust, vjust=vjust, label=sprintf("%i*",n)), size=3) +
        ggrepel::geom_label_repel(data = bind_rows(plot_data, missing),
            aes(label=label, color=type, fontface=fface),
            size=2, na.rm=TRUE, segment.alpha=0.3, fill="#ffffffc0",
            label.padding=0.1, max.iter=1e4, min.segment.length=0) +
        labs(subtitle = subt) +
        theme_classic()
}

sys$run({
    args = sys$cmd$parse(
        opt('d', 'dset', 'rds', 'pan.rds'),
        opt('c', 'cna', 'amp|del|all', 'all'),
        opt('f', 'fit', 'rank|rlm{,2,3}', 'rlm3'),
        opt('m', 'stat_max', 'numeric', '20'),
        opt('p', 'plotfile', 'pdf', 'cor_dset/pan_rlm2_amp.pdf'))

    assocs = readRDS(args$dset) %>%
        filter(dset == "orf" | (cna == args$cna & fit == args$fit)) %>%
        mutate(dset = sub("_none", "", paste(dset, adj, sep="_")))
    cap = as.numeric(args$stat_max)
    smat = narray::construct(statistic ~ name + dset, data=assocs)

    ds = unique(assocs$dset)
    plots = expand.grid(a1 = ds, a2 = ds, stringsAsFactors=FALSE) %>%
        filter(a1 < a2) %>%
        tbl_df() %>%
        mutate(plots = purrr::map2(a1, a2, do_plot, cap=cap, smat=smat,
                                   label=c(20, 10, 20, 3)))

    pdf(args$plotfile, 10, 10)
    for (i in seq_len(nrow(plots)))
        print(plots$plots[[i]] + ggtitle(sprintf("%s vs %s", plots$a1[i], plots$a2[i])))
    dev.off()
})
