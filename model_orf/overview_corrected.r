library(ggplot2)
io = import('io')
sys = import('sys')
ov = import('./overview_naive')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'overview.rds'),
    opt('p', 'plotfile', 'pdf', 'overview_corrected.pdf'))

expr = readRDS(args$infile)

percell = expr %>%
    mutate(cells = sprintf("%s (%s)", cells, tissue)) %>%
    group_by(cells) %>%
    tidyr::nest() %>%
    mutate(plot = purrr::map2(data, cells, ov$plot_overview_corrected))

pdf(args$plotfile)
for (p in percell$plot)
    print(p)
dev.off()
