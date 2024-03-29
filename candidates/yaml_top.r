library(dplyr)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'xlsx', '../merge_rank/rank_top/pan_rlm3.xlsx'),
    opt('n', 'num', 'integer', '12'),
    opt('o', 'outfile', 'yaml', 'pan/top-genes.yaml'))

n = as.integer(args$num)

sheets = readxl::excel_sheets(args$infile)
top = sapply(sheets, simplify=FALSE, function(s) {
    ranks = readxl::read_excel(args$infile, s)
})

res = list(
    genes = list(
        amp = top$amp$name[seq_len(n)],
        del = top$del$name[seq_len(n)]
    ),
    methods = c("lm", "stan-nb")
)
res$genes$both = setdiff(top$all$name, c(res$genes$amp,res$genes$del))[seq_len(n)]
yaml::write_yaml(res, file=args$outfile)
