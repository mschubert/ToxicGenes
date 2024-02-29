library(dplyr)
library(ggplot2)
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

both = inner_join(clin, rbm) |>
    mutate(gain = tcn_RBM14 >= 3)
colnames(both) = make.names(colnames(both))

m = coxph(Surv(Overall.survival.days, Vital.Status) ~ Sex + Age.at.diagnosis + Pre.Treated + Tumour.Grade + gain, data=both)
broom::tidy(m)

m2 = survfit(Surv(Overall.survival.days, Vital.Status) ~ Sex, data=both)
survplot(m2)
