library(dplyr)
library(ggplot2)
cm = import('../../report/common')

hct = readxl::read_xlsx("CRISPRa_Veronica/CRISPRAscreen_for Veronica_VR_HCT116_HahnLab.xlsx")

bt869 = list.files("CRISPRa_Veronica", pattern="LFC_DMSO_ETP.*", full.names=TRUE) |>
    lapply(readr::read_tsv) |>
    setNames(c("rep1", "rep2", "rep3")) |>
    bind_rows(.id="rep") |>
    group_by(gene = `Gene Symbol`) |>
    summarize(avg_lfc = mean(`Average LFC`))

argos = cm$get_argos(pan=TRUE)
tox = cm$get_tox()$`Pan-Cancer` |>
    inner_join(bt869) |>
    mutate(label = ifelse(gene %in% argos, gene, NA_character_))

m = broom::tidy(lm(avg_lfc ~ estimate, data=tox))

pdf("CRISPRa_Veronica.pdf")
ggplot(tox, aes(x=estimate, y=avg_lfc, color=is_toxic)) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label=label), color="black") +
    geom_smooth(method="lm", color="blue") +
    annotate("text", x=-1, y=1, color="blue",
             label=sprintf("P = %.2g", m$p.value[m$term != "(Intercept)"])) +
    ggtitle("BT869")
dev.off()
