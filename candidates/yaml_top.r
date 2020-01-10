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
        top12 = top$all$name[seq_len(n)],
        top24 = top$all$name[seq_len(2*n)[-seq_len(n)]],
        amp = top$amp$name[seq_len(n)],
        del = top$del$name[seq_len(n)]
    ),
    methods = c("lm", "rlm", "rlm3")
)
yaml::write_yaml(res, file=args$outfile)
