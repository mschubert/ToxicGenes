library(dplyr)
library(cowplot)
b = import('base')

dset = list(
        orf_pan = readxl::read_xlsx("../orf/fits_corrected.xlsx", "pan") %>%
            dplyr::rename(gene = `GENE SYMBOL`),
        orf_pancov = readxl::read_xlsx("../orf/fits_corrected.xlsx", "pancov") %>%
            dplyr::rename(gene = `GENE SYMBOL`),
        ccle = readxl::read_xlsx("../ccle/pan.xlsx", "amp"),
        tcga_naive = readxl::read_xlsx("../tcga/naive/pan.xlsx", "amp"),
        tcga_pur = readxl::read_xlsx("../tcga/pur/pan.xlsx", "amp"),
        tcga_puradj = readxl::read_xlsx("../tcga/pur+adj/pan.xlsx", "amp")
    )

assocs = dset %>%
    dplyr::bind_rows(.id="assocs") %>%
    select(-`Construct IDs`, -n_aneup) #%>%
#    filter(statistic < 0)

cap = 40

smat = assocs %>%
    group_by(gene, assocs) %>%
    arrange(-statistic) %>%
    top_n(1, "statistic") %>%
    ungroup() %>%
    mutate(statistic = sign(statistic) * pmin(abs(statistic), cap)) %>%
    dplyr::distinct(assocs, gene, .keep_all=TRUE) %>%
    narray::construct(statistic ~ gene + assocs)

do_plot = function(a1, a2) {
    x1 = rlang::sym(a1)
    x2 = rlang::sym(a2)
    colors = setNames(c("#e34a33", "#2ca25f", "#2b8cbe", "#cccccc"),
                      c("compensated", "hyper-dereg", "inconsistent", "no change"))
    wald = 1.5

    data = as.data.frame(na.omit(smat[,c(a1, a2)])) %>%
        tibble::rownames_to_column("gene")
    plot_data = data %>%
        mutate(type = case_when(
                abs(!! x1) < wald | abs(!! x2) < wald ~ "no change",
                !! x1 < 0 & !! x2 < 0 ~ "compensated",
                !! x1 > 0 & !! x2 > 0 ~ "hyper-dereg",
                TRUE ~ "inconsistent"),
               type = factor(type, levels=names(colors))) %>%
        group_by(type) %>%
        mutate(annot_score = abs(!! x1 / max(abs(!! x1), na.rm=TRUE) *
                                 !! x2 / max(abs(!! x2), na.rm=TRUE)),
               label = ifelse(rank(-annot_score) <= 20, gene, NA)) %>%
        ungroup()
    pcor = broom::tidy(lm(plot_data[[a1]] ~ plot_data[[a2]]))$p.value[2]

    nums = data %>%
        group_by(a1=sign(!! x1), a2=sign(!! x2)) %>%
        summarize(n = n()) %>%
        mutate(hjust=-a1/1.8+0.5, vjust=-a2+0.5)
    pfet = broom::tidy(fisher.test(matrix(pull(nums, n), ncol=2))) %catch%
        list(estimate=NA, p.value=NA)
    nums2 = data %>%
        group_by(a1=sign(!! x1), a2=sign(!! x2)) %>%
        filter(abs(!! x1) >= wald & abs(!! x2) >= wald) %>%
        summarize(n = n()) %>%
        mutate(x=sign(a1)*2, y=sign(a2)*2, hjust=-a1/2+0.5, vjust=-a2/2+0.5)
    pfet2 = broom::tidy(fisher.test(matrix(pull(nums2, n), ncol=2))) %catch%
        list(estimate=NA, p.value=NA)

    subt = sprintf(paste(sep = " | ",
                         "lm (red line): p=%.2g",
                         "comp. FET OR %.2f* (%.2f) p=%.2g* (%.2g)",
                         "stat cap %g",
                         "wald %g"),
                   pcor, pfet2$estimate, pfet$estimate,
                   pfet2$p.value, pfet$p.value, cap, wald)

    ggplot(plot_data, aes_string(x=a1, y=a2)) +
        geom_point(aes(color=type)) +
        geom_hline(yintercept = 0, size=1, linetype="dashed", color="#dedede") +
        geom_vline(xintercept = 0, size=1, linetype="dashed", color="#dedede") +
        scale_color_manual(values=colors) +
        geom_smooth(method="lm", color="#a50f15", linetype="dotted", se=FALSE) +
        geom_text(data=nums, aes(hjust=hjust, vjust=vjust, label=n), x=0, y=0, size=3) +
        geom_text(data=nums2, aes(x=x, y=y, hjust=hjust, vjust=vjust, label=sprintf("%i*",n)), size=3) +
        ggrepel::geom_label_repel(aes(label=label), size=2, na.rm=TRUE,
                                  fill="#ffffff70", label.padding=unit(0.1,"lines")) +
        labs(subtitle = subt)
}

plots = expand.grid(a1 = names(dset), a2 = names(dset), stringsAsFactors=FALSE) %>%
    filter(a1 < a2,
           ! a2 %in% c("tcga_naive", "tcga_pur"),
           a1 != "orf_pancov",
           ! (a1 == "ccle" & a2 == "orf_pancov")) %>%
    tbl_df() %>%
    mutate(plots = purrr::map2(a1, a2, do_plot))

pdf("cor.pdf", 10, 10)
for (i in seq_len(nrow(plots)))
    print(plots$plots[[i]] + ggtitle(sprintf("%s vs %s", plots$a1[i], plots$a2[i])))
dev.off()
