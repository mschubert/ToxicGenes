library(dplyr)

fra = list.files("./fragile_site_bed", full.names=TRUE) %>%
    lapply(readr::read_tsv, col_names=c("seqnames", "start", "end", "fra", "x", "strand")) %>%
    bind_rows() %>%
    select(-x, -strand) %>%
    mutate(seqnames = sub("chrx", "chrX", seqnames, fixed=TRUE)) %>%
    arrange(seqnames)

write.table(fra, file="fra.txt", sep="\t", quote=FALSE, row.names=FALSE)
