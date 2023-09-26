library(dplyr)
library(ggplot2)

orf_ccle = c(
#    H2077 = "",
#    BE2C = "",
#    NB69 = "",
#    `LAN-1` = "",
    Meljuso = "MELJUSO_SKIN",
    OVCAR4 = "OVCAR4_OVARY",
    OVSAHO = "OVSAHO_OVARY",
    Kura = "KURAMOCHI_OVARY", #??
    BT474 = "BT474_BREAST",
    D283 = "D283MED_CENTRAL_NERVOUS_SYSTEM",
#    D458 = "",
    T47D = "T47D_BREAST",
    SKBr3 = "SKBR3_BREAST",
    LnCaP = "LNCAPCLONEFGC_PROSTATE",
    TC32 = "TC32_BONE",
#    `SK-NEP-1` = "",
    `WM266-4` = "WM2664_SKIN"
)

ab = readxl::read_xlsx("../data/ccle/CCLE_ABSOLUTE_combined_20181227.xlsx") %>%
    filter(sample %in% orf_ccle) %>%
    group_by(sample) %>%
    summarize(ABS_ploidy = weighted.mean(Modal_Total_CN, Length),
              ABS_aneup = weighted.mean(abs(Modal_Total_CN - 2), Length))

dset = readRDS("overview.rds") %>%
    filter(`GENE SYMBOL` == "RBM14", `BEST GENE MATCH` == "RBM14") %>%
    mutate(sample = orf_ccle[cells]) %>%
    select(sample, tissue, ETP, `LFC DMSO/ETP`, z_LFC) %>%
    inner_join(ab)

long = dset %>%
    tidyr::pivot_longer(c(`LFC DMSO/ETP`, z_LFC), names_to=c("LFC_type"), values_to=c("LFC")) %>%
    tidyr::pivot_longer(c(ABS_ploidy, ABS_aneup), names_to=c("aneup_type"), values_to=c("aneup_value"))

ggplot(long, aes(x=aneup_value, y=LFC)) +
    geom_point(aes(color=tissue, size=ETP)) +
    geom_point(data=long %>% filter(is.na(ETP)), aes(color=tissue), size=1) +
    geom_text(data=long %>% filter(is.na(ETP)), label="ETP?", size=2) +
    geom_smooth(method="lm") +
    facet_grid(LFC_type ~ aneup_type, scales="free")
