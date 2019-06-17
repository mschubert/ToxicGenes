io = import('io')

expr = io$read_table("./data/ORF_DMSO_2019-02.txt", header=TRUE) %>%
    tidyr::gather("condition", "value", -(`Construct Barcode`:`BEST GENE MATCH`)) %>%
    mutate(condition = sub("( ORF)?[_-]DMSO", " DMSO", condition),
           condition = sub("LFC$", "LFC DMSO/ETP", condition),
           cells = sub(" .*$", "", condition),
           condition = sub("^[A-Za-z0-9-]+ ", "", condition)) %>%
    tidyr::spread(condition, value)
