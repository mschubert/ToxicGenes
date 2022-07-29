sys = import('sys')
idmap = import('process/idmap')
gset = import('genesets')

args = sys$cmd$parse(
    opt('h', 'haplo', 'xlsx' , 'pnas.1900437116.sd01.xlsx'),
    opt('e', 'essential', 'text file', 'depmap_essential.txt'),
    opt('n', 'nonessential', 'text file', 'depmap_nonessential.txt'),
    opt('d', 'davoli', 'xlsx', '1-s2.0-S0092867413012877-mmc2.xlsx'),
    opt('o', 'outfile', 'save to rds', 'manual.rds')
)

davoli = readxl::read_xlsx(args$davoli, skip=2)
cor = gset$get_human(c("CORUM_core", "CORUM_all", "CORUM_splice"))

cosmic = readr::read_tsv("Census_COSMIC96.tsv")
cosmic_og = cosmic %>% filter(grepl("oncogene", `Role in Cancer`))
cosmic_tsg = cosmic %>% filter(grepl("TSG", `Role in Cancer`))

sets = list(
    haplo_insuff = setdiff(readxl::read_xlsx(args$haplo, "published data")$gene, "WT"),
    DepMap_essential = sub("^([^ ]+).*", "\\1",
            read.table(args$essential, header=TRUE, sep="\t")$gene),
    DepMap_nonessential = sub("^([^ ]+).*", "\\1",
            read.table(args$nonessential, header=TRUE, sep="\t")$gene),
    Davoli_oncogenes = davoli$OG[!is.na(davoli$OG)],
    Davoli_TSGs = davoli$TSG[!is.na(davoli$TSG)],
    Cosmic_OG_Hallmark = cosmic_og %>% filter(Tier == 1, Hallmark == "Yes") %>% pull(`Gene Symbol`),
    Cosmic_OG_Hallmark = cosmic_tsg %>% filter(Tier == 1, Hallmark == "Yes") %>% pull(`Gene Symbol`),
    Cosmic_OG_Tier1 = cosmic_og %>% filter(Tier == 1) %>% pull(`Gene Symbol`),
    Cosmic_OG_Tier2 = cosmic_og %>% filter(Tier == 2) %>% pull(`Gene Symbol`),
    Cosmic_TSG_Tier1 = cosmic_tsg %>% filter(Tier == 1) %>% pull(`Gene Symbol`),
    Cosmic_TSG_Tier2 = cosmic_tsg %>% filter(Tier == 2) %>% pull(`Gene Symbol`),
    CORUM_all= unique(unlist(cor$CORUM_all)),
    CORUM_core= unique(unlist(cor$CORUM_core)),
    CORUM_splice= unique(unlist(cor$CORUM_splice))
)

sets = c(
    sets,
    cor$CORUM_all[grepl("snRNP|proteasome", names(cor$CORUM_all)) |
                  sapply(cor$CORUM_all, function(x) "RBM14" %in% x)]
)

saveRDS(sets, file=args$outfile)
