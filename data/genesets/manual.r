sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('h', 'haplo', 'xlsx' , 'pnas.1900437116.sd01.xlsx'),
    opt('e', 'essential', 'text file', 'depmap_essential.txt'),
    opt('n', 'nonessential', 'text file', 'depmap_nonessential.txt'),
    opt('d', 'davoli', 'xlsx', '1-s2.0-S0092867413012877-mmc2.xlsx'),
    opt('o', 'outfile', 'save to rds', 'manual.rds'))

davoli = readxl::read_xlsx(args$davoli, skip=2)

sets = list(
    haplo_insuff = setdiff(readxl::read_xlsx(args$haplo, "published data")$gene, "WT"),
    common_essential = sub("^([^ ]+).*", "\\1",
            read.table(args$essential, header=TRUE, sep="\t")$gene),
    non_essential = sub("^([^ ]+).*", "\\1",
            read.table(args$nonessential, header=TRUE, sep="\t")$gene),
    oncogenes = davoli$OG,
    TSGs = davoli$TSG
)

saveRDS(sets, file=args$outfile)