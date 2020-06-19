io = import('io')
sys = import('sys')

score = function(dset) {
    score = dset %>%
        filter(adj %in% c("none", "pur"), #TODO: keep "none" for TCGA too???
               fit %in% c("lm", args$fit)) %>%
        mutate(rsq = ifelse(dset == "orf", 0.5, rsq), #TODO: better way? (value determines ORF weight)
               estimate = ifelse(dset == "orf", 2^estimate - 1, estimate),
               estimate = pmax(estimate, -1)) %>%
        mutate(score = (1-adj.p) * (rank(-statistic) / length(statistic)) *
                       pmax(0,rsq) * (-estimate))

    score %>%
        group_by(name, dset) %>% # mean by fit (rlm, rank)
        summarize(score = mean(score, na.rm=TRUE)) %>%
        group_by(name) %>% # summarize by dset (orf, ccle, tcga)
        summarize(n_dset = length(name),
                  score = sum(score^0.5, na.rm=TRUE)) %>% # root prioritizes consistency across dsets
        arrange(-score) %>%
        mutate(rank = 1:nrow(.))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'dset', '../merge/pan.rds'),
    opt('f', 'fit', 'character', 'rlm3'),
    opt('o', 'outfile', 'xlsx', 'pan.xlsx'))

dset = readRDS(args$infile)
top = list(
    all = score(dset),
    amp = score(dset %>% filter(cna %in% c("oe", "amp"))),
    del = score(dset %>% filter(cna == "del"))
)
writexl::write_xlsx(top, args$outfile)
