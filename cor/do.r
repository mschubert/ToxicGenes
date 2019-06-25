library(dplyr)
library(cowplot)

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

amat = assocs %>%
    group_by(gene, assocs) %>%
    arrange(-statistic) %>%
    top_n(1, "statistic") %>%
    ungroup() %>%
    mutate(statistic = sign(statistic) * pmin(abs(statistic), cap)) %>%
    dplyr::distinct(assocs, gene, .keep_all=TRUE) %>%
#    mutate(statistic = abs(statistic)) %>%
    narray::construct(statistic ~ gene + assocs)

do_plot = function(a1, a2) {
    x1 = rlang::sym(a1)
    x2 = rlang::sym(a2)
    data = as.data.frame(na.omit(amat[,c(a1, a2)])) %>%
        tibble::rownames_to_column("gene") %>%
        mutate(type = case_when(
                !! x1 < 0 & !! x2 < 0 ~ "compensated",
                !! x1 > 0 & !! x2 > 0 ~ "hyper-dereg",
                TRUE ~ "inconsistent")) %>%
        group_by(type) %>%
        mutate(annot_score = abs(!! x1 / max(abs(!! x1), na.rm=TRUE) *
                                 !! x2 / max(abs(!! x2), na.rm=TRUE)),
               label = ifelse(rank(-annot_score) <= 20, gene, NA)) %>%
        ungroup()
    pcor = broom::tidy(lm(data[[a1]] ~ data[[a2]]))$p.value[2]

    ggplot(data, aes_string(x=a1, y=a2)) +
        geom_hline(yintercept = 0, size=1, linetype="dashed", color="grey") +
        geom_vline(xintercept = 0, size=1, linetype="dashed", color="grey") +
        geom_point(aes(color=type)) +
        geom_smooth(method="lm", color="red", linetype="dotted", se=FALSE) +
        ggrepel::geom_label_repel(aes(label=label), size=2) +
        labs(subtitle = sprintf("correlation (red line): p=%.2g (stat capped at %i)", pcor, cap))
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
