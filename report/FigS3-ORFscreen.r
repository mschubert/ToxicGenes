library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')
orf = import('../orf/overview_naive')

facet_plot = function(ov, aes) {
    loess_sd = ov %>% group_by(cells) %>%
        summarize(mod = list(msir::loess.sd(!! aes$x, !! aes$y))) %>%
        rowwise() %>%
        mutate(res = list(data.frame(x=mod$x, y=mod$y, sd=mod$sd))) %>%
        select(-mod) %>%
        tidyr::unnest(res)

    ggplot(ov, aes) +
        geom_hline(yintercept=0, color="black") +
        geom_hex(aes(color=..count..), bins=50) +
        scale_color_continuous(type = "viridis", trans="log1p", guide="none") +
        scale_fill_continuous(type = "viridis", trans="log1p", breaks=c(1,5,20,100,500)) +
        geom_line(data=loess_sd, aes(x=x, y=y), color="red") +
        geom_line(data=loess_sd, aes(x=x, y=y+sd), color="red", linetype="dashed") +
        geom_line(data=loess_sd, aes(x=x, y=y-sd), color="red", linetype="dashed") +
        facet_wrap(~ cells, ncol=6) +
        theme_minimal() +
        theme(strip.background = element_rect(color=NA, fill="#ffffffc0"))
}

screen_cor = function(ov) {
    wide = ov %>%
        select(`Construct IDs`, cells, `LFC DMSO/ETP`) %>%
        tidyr::pivot_wider(names_from=cells, values_from=`LFC DMSO/ETP`)
    mat = data.matrix(wide[-1])
    rownames(mat) = wide[[1]]

    cmat = cor(mat) %>%
        reshape2::melt() %>%
        plt$cluster(value ~ Var1 + Var2)
    plt$matrix(cmat, value ~ Var1 + Var2, geom="tile") +
        scale_fill_distiller(palette="RdBu", name="Pearson\ncorrelation") +
        theme(axis.title = element_blank()) +
        coord_fixed()
}

go_volc = function() {
    res = readxl::read_xlsx("../orf/pan/GO_Biological_Process_2018.xlsx")
    plt$volcano(res, label_top=35) + guides(size="none") +
        xlab("Mean z-score LFC")
}

sys$run({
    ov = readRDS("../orf/overview.rds") %>%
        mutate(cells = sprintf("%s (%s)", cells, tissue))

    naive = facet_plot(ov, aes(x=DMSO, y=`LFC DMSO/ETP`)) +
        ylim(c(-4,4)) +
        coord_cartesian(ylim=c(-2.6,2.6), clip="off")
    corr = facet_plot(ov, aes(x=DMSO, y=z_LFC)) +
        ylim(c(-9,9)) +
        coord_cartesian(ylim=c(-6,6), clip="off")

    cors = wrap_elements(screen_cor(ov))

    asm = (((naive / corr) | ((cors / go_volc()) + plot_layout(heights=c(1,2)))) +
        plot_layout(widths=c(1.5,1))) + plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS3-ORFscreen.pdf", 14, 10)
    print(asm)
    dev.off()
})