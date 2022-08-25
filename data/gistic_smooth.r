library(dplyr)
sys = import('sys')
seq = import('seq')

get_gistic_scores = function(gistic) {
    gt = seq$gene_table() %>%
        transmute(chr = factor(chromosome_name, levels=c(1:22,'X')),
                  gene_name = external_gene_name,
                  tss = transcription_start_site) %>%
        filter(!is.na(chr)) %>%
        group_by(chr, gene_name) %>%
            summarize(tss = mean(tss)) %>%
        ungroup()

     gistic %>%
        transmute(gene_name = SYMBOL,
                  type = ALT_TYPE,
                  frac = ifelse(type == "amplification", OVERALL_FREQ, -OVERALL_FREQ)) %>%
        inner_join(gt)
}

cna2smooth = function(chr, steps) {
    fracs = chr %>% arrange(tss)
    gam = mgcv::gam(frac ~ s(tss), data=fracs)

    res = dplyr::tibble(tss = steps)
    res$freq = predict(gam, newdata=res)
    res
}

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'tcga_prod2-gistic.rds'),
    opt('o', 'outfile', 'rds', 'gistic_smooth.rds')
)

gistic = readRDS(args$infile)
scores = get_gistic_scores(gistic)
pred_x = seq$chr_step("GRCh38", step=1e6, chrs=c(1:22,'X'))

res = scores %>%
    group_by(type, chr) %>%
    tidyr::nest() %>%
    arrange(type, chr) %>%
    inner_join(pred_x) %>%
    rowwise() %>%
    mutate(gam = list(cna2smooth(data, steps)))

smooth = res %>%
    select(type, chr, gam) %>%
    tidyr::unnest(gam)

saveRDS(list(genes=scores, smooth=smooth), file=args$outfile)
