library(dplyr)

con = DBI::dbConnect(RMariaDB::MariaDB(), dbname = "tcga")
tabs = DBI::dbListTables(con) %>%
    sapply(function(x) tbl(con, x), simplify=FALSE)

tabs = tabs[sapply(tabs, function(x) tally(x) %>% pull(n) != 0)]
# ALT: alteration
# TSU: tissue

gstc_analysis = tbl(con, "GSTC_ANALYSIS") %>%
    filter(VERSION == "2015-06-01 stddata__2015_04_02 arm-level peel-off") %>%
    select(GSTC_ANALYSIS_ID, GSTC_METADATA_ID)
gstc_tsu = tbl(con, "GSTC_TSU") %>%
    filter(TISSUE_TYPE == "all_cancers") %>%
    select(GSTC_TSU_ID, GSTC_METADATA_ID)

gstc_gene_alt = tbl(con, "GSTC_GENE_ALT") %>%
    select(GSTC_GENE_ANALYSIS_ID, GSTC_GENE_ALT_ID, ALT_TYPE)
gstc_genome = tbl(con, "GSTC_GENOME") %>%
    filter(BUILD == "hg19") %>%
    select(GSTC_GENOME_ID)
gstc_gene = tbl(con, "GSTC_GENE") %>%
    inner_join(gstc_genome) %>%
    select(GSTC_GENE_ID, SYMBOL, REGION)

gene = tbl(con, "GSTC_GENE_ANALYSIS") %>%
    inner_join(gstc_gene) %>%
    inner_join(gstc_gene_alt)

meta = gstc_analysis %>%
    inner_join(gstc_tsu)

res = tbl(con, "GSTC_GENE_ALT_ANALYSIS") %>%
    inner_join(meta) %>%
    inner_join(gene) %>%
    select(SYMBOL, FOCAL_FREQ, HIGH_LEVEL_FREQ, OVERALL_FREQ, ALT_TYPE)

saveRDS(as_tibble(res), file="tcga_prod2-gistic.rds")
