library(dplyr)

total_quant = function(fname) {
    x = readr::read_tsv(fname) %>% arrange(FDR, PValue)

    inc_luc = strsplit(x$IJC_SAMPLE_1, ",") %>% do.call(rbind, .)
    skip_luc = strsplit(x$SJC_SAMPLE_1, ",") %>% do.call(rbind, .)
    inc_rbm = strsplit(x$IJC_SAMPLE_2, ",") %>% do.call(rbind, .)
    skip_rbm = strsplit(x$SJC_SAMPLE_2, ",") %>% do.call(rbind, .)

    mode(inc_luc) = "numeric"
    mode(skip_luc) = "numeric"
    mode(inc_rbm) = "numeric"
    mode(skip_rbm) = "numeric"

    luc = colSums(skip_luc) / (colSums(inc_luc / 2) + colSums(skip_luc))
    rbm = colSums(skip_rbm) / (colSums(inc_rbm / 2) + colSums(skip_rbm))
    rbind(luc=luc, rbm=rbm)
}

total_quant("rmats_out/RI.MATS.JC.txt")
total_quant("rmats_out/MXE.MATS.JC.txt")
total_quant("rmats_out/A3SS.MATS.JC.txt")
total_quant("rmats_out/A5SS.MATS.JC.txt")
