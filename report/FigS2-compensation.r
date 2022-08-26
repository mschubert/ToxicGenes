library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')

sys$run({
    asm = plot_spacer()

    asm = norms + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS2-compensation.pdf", 14, 16)
    print(asm)
    dev.off()
})
