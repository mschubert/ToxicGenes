library(dplyr)
library(ggplot2)
library(patchwork)
library(survival)

clin = readxl::read_xlsx("Supplementary_Table_01.xlsx", skip=2) |>
    mutate(sample = `DNA Tumor Sample Barcode`,
           `Age at diagnosis` = as.numeric(`Age at diagnosis`),
           `Vital Status` = as.integer(factor(`Vital Status`, levels=c("Alive", "Dead")))-1,
           `Overall survival days` = as.numeric(`Overall survival days`)) |>
    filter(!is.na(`Overall survival days`))
cns = readr::read_tsv("CRC-SW.FACETS.facets-suite.1063_DNBSEQ.20210706.gene_level.gencode_v35.ProteinCoding_IG_TR.txt")

rbm = cns |>
    filter(gene %in% c("CCND1", "RBM14")) |>
    select(sample, gene, median_cnlr_seg, tcn) |>
    tidyr::pivot_wider(names_from="gene", values_from=c("median_cnlr_seg", "tcn"))

rna = readr::read_tsv("CRC.SW.mRNA.symbol.TPM.txt.gz") |>
    filter(SYMBOL %in% c("CCND1", "RBM14")) |>
    tidyr::pivot_longer(-SYMBOL, names_to="sample", values_to="tpm") |>
    tidyr::pivot_wider(names_from=SYMBOL, values_from="tpm", names_prefix="tpm_") |>
    mutate(sample = gsub(".", "-", sample, fixed=TRUE))

both = inner_join(clin, rbm) |>
    left_join(rna) |>
    mutate(
        gainCCND1 = case_when(
            tcn_CCND1 >= 4 ~ TRUE,
            tcn_CCND1 == 2 ~ FALSE,
            TRUE ~ NA
        ),
        gainRBM14 = case_when(
            tcn_RBM14 >= 4 ~ TRUE,
            tcn_RBM14 == 2 ~ FALSE,
            TRUE ~ NA
        ),
        group = case_when(
            !gainCCND1 & gainRBM14 ~ NA_character_, # 1 sample
            gainCCND1 & gainRBM14 ~ "CCND1_RBM14",
            gainCCND1 ~ "CCND1",
            TRUE ~ "neither"
        ),
        group = factor(group, levels=c("CCND1", "CCND1_RBM14", "neither"))
    )
colnames(both) = make.names(colnames(both))

m = coxph(Surv(Overall.survival.days, Vital.Status) ~ Sex + Age.at.diagnosis + group, data=both)
broom::tidy(m)

m2 = survfit(Surv(Overall.survival.days, Vital.Status) ~ group, data=both)
m2
p = ggsurvplot(m2, data=both)

p2 = ggplot(both, aes(x=tcn_CCND1, y=tpm_CCND1)) +
    geom_point() +
    geom_smooth(method="lm") +
    xlim(0,10) + ylim(0, 650)

p3 = ggplot(both, aes(x=tcn_RBM14, y=tpm_RBM14)) +
    geom_point() +
    geom_smooth(method="lm") +
    xlim(0,NA) + ylim(0,NA)

pdf()
print(p)
print(p2 | p3)
dev.off()
