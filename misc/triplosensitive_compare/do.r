library(dplyr)
library(betareg)
library(ggplot2)
library(patchwork)

comp = readr::read_tsv("../../cor_tcga_ccle/positive_comp_set.tsv")
compg = comp$gene[comp$hit]
argos =  comp$gene[comp$hit & !is.na(comp$stat_orf) & comp$stat_orf < -5]

dset = readxl::read_xlsx("1-s2.0-S0092867422007887-mmc7.xlsx") #%>%
#    filter(pTriplo > 0.9)

dset$is_comp = dset$Gene %in% compg
dset$is_argos = dset$Gene %in% argos

broom::tidy(betareg(pTriplo ~ is_comp, data=dset))

p1 = ggplot(dset, aes(x=is_comp, y=pTriplo)) +
    ggbeeswarm::geom_quasirandom(alpha=0.5) +
    ggpubr::stat_compare_means() +
    stat_summary(fun=mean, geom="crossbar", colour="red")

p2 = ggplot(dset, aes(x=is_argos, y=pTriplo)) +
    ggbeeswarm::geom_quasirandom(alpha=0.5) +
    ggpubr::stat_compare_means() +
    stat_summary(fun=mean, geom="crossbar", colour="red")

(p1 | p2) & ylim(NA, 1.05)
