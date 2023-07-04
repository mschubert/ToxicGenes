library(dplyr)
library(ggplot2)
library(patchwork)
library(survival)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'plot_cnvs.rds'),
    opt('m', 'meta', 'xlsx', 'N4plus_data_Michael_20230414.xlsx'),
    opt('p', 'plotfile', 'pdf', 'surv.pdf')
)

meta = readxl::read_xlsx(args$meta) %>%
    mutate_all(function(x) ifelse(x == 9999, NA, x)) %>%
    mutate(PtID = sprintf("%03i", PtID))

dset = readRDS(args$infile)
segs = dset$segs %>%
    filter(chr == "chr11") %>%
    inner_join(meta %>% select(PtID, TCP)) %>%
    mutate(mean = log2(1+(2^mean-1)/(TCP/100)))

genes = tibble(genes = c("CCND1", "RBM14"), pos=c(69647000, 66622000))

p1 = ggplot(segs) +
    geom_hline(yintercept=0, size=1, color="red") +
    geom_segment(aes(x=start, xend=end, y=mean, yend=mean), alpha=0.1) +
    geom_point(data=genes, aes(x=pos), y=0, color="red", size=3, alpha=0.5) +
    ggrepel::geom_label_repel(data=genes, aes(x=pos, label=genes), y=0,
                              color="red", alpha=0.5) +
    coord_cartesian(expand=FALSE)

segs2 = segs %>% filter(start < 69647000, end > 66622000)
segs3 = segs2 %>% filter(start < 66622000, end > 66622000) # RBM14
segs4 = segs2 %>% filter(start < 69647000, end > 69647000) # CCND1
p2 = ggplot(segs2) +
    geom_boxplot(data=segs3, aes(y=mean), x=66622000, fill="white", alpha=0.5) +
    geom_boxplot(data=segs4, aes(y=mean), x=69647000, fill="white", alpha=0.5) +
    geom_hline(yintercept=0, size=1, color="red") +
    geom_segment(aes(x=start, xend=end, y=mean, yend=mean), alpha=0.1) +
    geom_point(data=genes, aes(x=pos), y=0, color="red", size=3, alpha=0.5) +
    ggrepel::geom_label_repel(data=genes, aes(x=pos, label=genes), y=0,
                              color="red", alpha=0.5) +
    coord_cartesian(expand=FALSE, xlim=c(5e7,9e7))

meta2 = meta %>%
    mutate(CCND1 = case_when(PtID %in% segs4$PtID[segs4$mean>0.5] ~ "amp",
                             PtID %in% segs4$PtID[abs(segs4$mean)<0.1] ~ "eu"),
           RBM14 = case_when(PtID %in% segs3$PtID[segs3$mean>0.5] ~ "amp",
                             PtID %in% segs3$PtID[abs(segs3$mean)<0.1] ~ "eu"))
s1 = survfit(Surv(osDays_FINAL, OS_event_FINAL) ~ RBM14 + Treatment, data=meta2)
s2 = survfit(Surv(osDays_FINAL, OS_event_FINAL) ~ CCND1 + Treatment, data=meta2)
m11 = broom::tidy(coxph(Surv(osDays_FINAL, OS_event_FINAL) ~ T_ER + T_LN10 + Treatment,
                  data=meta2 %>% filter(RBM14 == "amp")))
m12 = broom::tidy(coxph(Surv(osDays_FINAL, OS_event_FINAL) ~ T_ER + T_LN10 + Treatment,
                  data=meta2 %>% filter(RBM14 == "eu")))
m21 = broom::tidy(coxph(Surv(osDays_FINAL, OS_event_FINAL) ~ T_ER + T_LN10 + Treatment,
                  data=meta2 %>% filter(CCND1 == "amp")))
m22 = broom::tidy(coxph(Surv(osDays_FINAL, OS_event_FINAL) ~ T_ER + T_LN10 + Treatment,
                  data=meta2 %>% filter(CCND1 == "eu")))

make_title = function(main, m1, m2) {
    m1 = m1 %>% filter(term == "Treatment")
    m2 = m2 %>% filter(term == "Treatment")
    sprintf("%s p=%2.g (amp), p=%.2g (eu)", main, m1$p.value, m2$p.value)
}

pdf(args$plotfile, 9, 6)
print(p1)
print(p2)
ggsurvplot(s1, title=make_title("RBM14", m11, m12))
ggsurvplot(s2, title=make_title("CCND1", m21, m22))
dev.off()
