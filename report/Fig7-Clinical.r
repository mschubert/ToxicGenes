library(dplyr)
library(ggplot2)
library(patchwork)
library(survival)
sys = import('sys')
seq = import('seq')
plt = import('plot')
cm = import('./common')

sys$run({
    both = readRDS("../misc/crc_translation/do.rds")

    m = coxph(Surv(Overall.survival.days, Vital.Status) ~ Sex + Age.at.diagnosis + group, data=both)
    broom::tidy(m)

    m2 = survfit(Surv(Overall.survival.days, Vital.Status) ~ group, data=both)
    m2
    p = survminer::ggsurvplot(m2, data=both)

    cols = c(neither="#00bfc4", CCND1="#f8766d", CCND1_RBM14="#7cae00")

    p2 = ggplot(both, aes(x=tcn_CCND1, y=tpm_CCND1)) +
        geom_point(aes(fill=group), position=position_jitter(width=0.2), shape=21, alpha=0.8) +
        scale_fill_manual(values=cols, na.translate=FALSE) +
        geom_smooth(method="lm") +
        xlim(0,6) + ylim(0, 650) +
        labs(x = "Copy number CCND1",
             y = "Gene expression CCND1 (tpm)") +
        cm$theme_minimal()

    p3 = ggplot(both, aes(x=tcn_RBM14, y=tpm_RBM14)) +
        geom_point(aes(fill=group), position=position_jitter(width=0.2), shape=21, alpha=0.8) +
        scale_fill_manual(values=cols, na.translate=FALSE) +
        geom_smooth(method="lm") +
        xlim(0,6) + ylim(0,180) +
        labs(x = "Copy number RBM14",
             y = "Gene expression RBM14 (tpm)") +
        cm$theme_minimal()

    asm = ((plot_spacer() + p2 + p3 + (p$plot + guides(color=FALSE))) & cm$text_sizes()) +
        plot_annotation(tag_levels='a') + plot_layout(ncol=2, nrow=2, guides="collect")

    pdf("Fig7-Clinical.pdf", 10, 8)
    print(asm)
    dev.off()
})
